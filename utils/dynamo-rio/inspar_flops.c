/* Instruction Parallelism:
 * inspar.c
 *
 * The approach taken in this implementation does not compute the
 * dependency graph. Instead, dependencies are calculated using a much
 * simpler algorithm which just considers the dependencies for each
 * register. For data dependency analysis, we actually only need
 * to know when each register was last written to and read from.
 *
 * In the dependency graph, the different sets of parallel schedulable
 * instructions are separated by determining the frontiers of the graph.
 * However, we can also calculate them directly. The set an instruction
 * belongs to can be determined by calculating the sets for each of the
 * hazards, and then using the maximum of these values to ensure no hazards
 * will be present in the parallel version of the code.
 *
 * For each register referenced in an instruction, we calculate the set number
 * based on the three simple rules:
 *
 *   - RaW: The set # of an instruction that READS register X must be at
 *          least 1 greater than the highest WRITTEN TO set # for register X
 *
 *   - WaR: The set # of an instruction that WRITES register X must be at
 *          least 1 greater than the highest READ FROM set # for register X
 *
 *   - WaW: The set # of an instruction that WRITES register X must be at
 *          least 1 greater than the highest WRITTEN TO set # for register X
 *
 *   - (RaR produces no rule.)
 *
 * Now, the WRITTEN TO and READ FROM values of all registers referenced in the
 * instruction are updated to the set # of the current instruction.
 *
 * Lets create a table to visualize these rules. The table contains the calculated
 * set numbers for READ FROM and WRITTEN TO for each register.
 *
 *                           reg1   reg2   reg3   reg4   reg5   eflags
 *   instruction             r/w    r/w    r/w    r/w    r/w     r/w
 *   ===========           ============================================
 *   a. reg1 = 50              w1 |      |      |      |      |         -> set 1 = {a}
 *   b. reg2 = 20                 |   w1 |      |      |      |         -> set 1 = {a,b}
 *   c. reg3 = reg1 + 1     r2    |      |   w2 |      |      |         -> set 2 = {c}
 *   d. reg4 = reg1 + reg2  r2    |r2    |      |   w2 |      |         -> set 2 = {c,d}
 *   e. reg2 = reg4               |   w3 |      |r3    |      |         -> set 3 = {e}
 *   f. reg4 = 10                 |      |      |   w4 |      |         -> set 4 = {f}
 *   g. reg2 = reg1         r4    |   w4 |      |      |      |         -> set 4 = {f,g}
 *   h. reg5 = reg1         r2!   |      |      |      |   w2 |         -> set 2 = {c,d,h}
 *   i. test(reg3,reg1)     r3!   |      |r3    |      |      |    w3   -> set 3 = {e,i}
 *   j. jnz                       |      |      |      |      |r4       -> set 4 = {f,g,j}
 *
 *  (Notice the special case that is visible for reg1 in instruction h: The
 *  calculated value is 2 for READ FROM reg1, but the table already contained
 *  a higher value of 4 for READ FROM reg1. In the actual implementation, only
 *  the maximum value is stored, eliminating the need to search for the highest
 *  values)
 *
 *  Now that we have all the parallel schedulable sets, we can calculate the ILP:
 *   set   instructions size(set)
 *   ===   ============ ========
 *    1     {a,b}          2     \
 *    2     {c,d,h}        3      \
 *    3     {e,i}          2       } ILP = (2+3+2+3)/4 = 10/4 = 2.5
 *    4     {f,g,j}        3      /
 *
 *  But not even this is necessary! As the total number of instructions is always
 *  the same, ILP can just be calculated as:
 *
 *    ILP = total instructions / number of sets
 *
 * Author: Jeremie Raymond Miserez, National University of Singapore
 *        <a01068654@nus.edu.sg>, <jeremie@miserez.org>
 * */

#include <stddef.h> /* for offsetof */
#include "dr_api.h"
#include "dr_ir_instr.h"
#include <float.h> /* for FLT_MAX */

#ifdef WINDOW
# define DISPLAY_STRING(msg) dr_messagebox(msg)
#else
# define DISPLAY_STRING(msg) dr_printf("%s\n", msg);
#endif

#define NULL_TERMINATE(buf) buf[(sizeof(buf)/sizeof(buf[0])) - 1] = '\0'

#define TESTALL(mask, var) (((mask) & (var)) == (mask))
#define TESTANY(mask, var) (((mask) & (var)) != 0)

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAX3(a,b,c) (MAX(MAX((a),(b)),(c)))
#define MAX4(a,b,c,d) (MAX(MAX(MAX((a),(b)),(c)),(d)))

#define MY_MAX_BB 1000000
#define MY_NUM_EFLAGS 11
#define MY_EFLAGS_OFFSET DR_REG_LAST_VALID_ENUM

#define CONSIDER_MORE_REGS
#define CONSIDER_EFLAGS
//#define USE_CLEAN_CALL
//#define INSERT_AT_END // does not work properly with firefox and multithreaded programs

/**
  */

struct srct_glob_reg_state {
    int raw_setnr;
    int war_setnr;
    int waw_setnr;
    int num_sets;
    int final_setnr;
    int else_setnr;
    int *my_readfrom;
    int *my_writtento;
};

typedef struct srct_glob_reg_state t_glob_reg_state;

/* Globals */
static void *stats_mutex; /* for multithread support */
static int my_bbcount;
static int my_bbexecs[MY_MAX_BB];
static int my_bbsizes[MY_MAX_BB];
static float my_bbilp[MY_MAX_BB];
static int bb_flop_count[MY_MAX_BB];

static void event_exit(void);
static dr_emit_flags_t event_basic_block(void *drcontext, void *tag, instrlist_t *bb,
                                         bool for_trace, bool translating);
inline void update_setnrs(instr_t *instr, int dr_reg_enum_start, int dr_reg_enum_stop,
                          t_glob_reg_state* glob_reg_state);
inline void update_eflag_setnrs(instr_t *instr, t_glob_reg_state* glob_reg_state);
#ifdef USE_CLEAN_CALL
static void clean_call(uint cur_num);
#endif

DR_EXPORT void 
dr_init(client_id_t id)
{
    dr_printf("dr_init()\n");
    my_bbcount = 0;
    stats_mutex = dr_mutex_create();
    dr_register_exit_event(event_exit);
    dr_register_bb_event(event_basic_block);

#ifdef SHOW_RESULTS
    dr_log(NULL, LOG_ALL, 1, "Client 'inspar' initializing\n");
    if (dr_is_notify_on()) {
# ifdef WINDOWS
        dr_enable_console_printing();
# endif
        dr_fprintf(STDERR, "Client inspar is running\n");
    }
#endif
}

static void 
event_exit(void)
{
    unsigned long long flop_count = 0;
    double ilp;
    int total_bbexecs = 0;
    double total_ilp = 0;
    float max_ilp = FLT_MIN;
    float min_ilp = FLT_MAX;
    int i = 0;
    for(i = 0; i < MY_MAX_BB && i < my_bbcount; i++){
        flop_count += my_bbexecs[i] * bb_flop_count[i];
        total_bbexecs += my_bbexecs[i];
        total_ilp += ((double)my_bbexecs[i])*((double)my_bbilp[i]);
        min_ilp = MIN(min_ilp, my_bbilp[i]);
        max_ilp = MAX(max_ilp, my_bbilp[i]);
    }
    ilp = total_ilp / ((double)(total_bbexecs != 0 ? total_bbexecs : 1));

#ifdef SHOW_RESULTS
    char msg[2048];
    int len;
#ifdef VERBOSE
    for(i = 0; i < MY_MAX_BB && i < my_bbcount; i++){
        len = dr_snprintf(msg, sizeof(msg)/sizeof(msg[0]),
                          "Basic block # %5d, size: %5d, execs: %8d, ILP: %3f",
                          i, my_bbsizes[i], my_bbexecs[i], my_bbilp[i]);
        DR_ASSERT(len > 0);
        NULL_TERMINATE(msg);
        dr_printf("%s\n", msg);
    }
#endif
    len = dr_snprintf(msg, sizeof(msg)/sizeof(msg[0]),
                      "Instrumentation results:\n"
		      "%10d basic block execs\n"
                      "%10d basic blocks\n"
            "%10d flops\n"
                      "%10.3f avg ILP\n"
                      "%10.3f max ILP\n"
                      "%10.3f min ILP\n",
                      total_bbexecs, my_bbcount, flop_count, ilp, max_ilp, min_ilp);
    DR_ASSERT(len > 0);
    NULL_TERMINATE(msg);
    dr_printf("%s\n", msg);
    if (my_bbcount > MY_MAX_BB){
        len = dr_snprintf(msg, sizeof(msg)/sizeof(msg[0]),
                          "\n"
                          "Overflow warning! Only the latest %d basic blocks were used for ILP calcution.\n",
                          MY_MAX_BB);
        DR_ASSERT(len > 0);
        NULL_TERMINATE(msg);
        dr_printf("%s\n", msg);
    }

#endif /* SHOW_RESULTS */
    dr_mutex_destroy(stats_mutex);
}


static inline void
calc_set_num(instr_t *instr, t_glob_reg_state* glob_reg_state )
{
  int i;
  //general purpose registers. Instead of using the whole enum, we only need to check for the
  //largest registers, any overlapping accesses will also trigger.
#ifdef X64
  update_setnrs(instr, DR_REG_START_64, DR_REG_STOP_64, glob_reg_state);
#else
  update_setnrs(instr, DR_REG_START_32, DR_REG_STOP_32, glob_reg_state);
#endif
#ifdef CONSIDER_MORE_REGS
  update_setnrs(instr, DR_REG_START_MMX, DR_REG_STOP_MMX, glob_reg_state);
  update_setnrs(instr, DR_REG_START_XMM, DR_REG_STOP_XMM, glob_reg_state);
  update_setnrs(instr, DR_REG_START_YMM, DR_REG_STOP_YMM, glob_reg_state);
  update_setnrs(instr, DR_REG_START_FLOAT, DR_REG_STOP_FLOAT, glob_reg_state);
  update_setnrs(instr, DR_REG_START_SEGMENT, DR_REG_START_CR, glob_reg_state);
  update_setnrs(instr, DR_REG_START_DR, DR_REG_STOP_DR, glob_reg_state);
  update_setnrs(instr, DR_REG_START_CR, DR_REG_STOP_CR, glob_reg_state);
#endif
#ifdef CONSIDER_EFLAGS
  update_eflag_setnrs(instr, glob_reg_state);
#endif
  glob_reg_state->final_setnr = MAX4(glob_reg_state->raw_setnr, glob_reg_state->war_setnr,
                        glob_reg_state->waw_setnr, glob_reg_state->else_setnr);
  //assigning the set number
  for (i = 0; i <= DR_REG_LAST_VALID_ENUM; i++) {
    if ( (i != DR_REG_NULL) && (i != DR_REG_INVALID) ) {
        if(instr_reads_from_reg(instr, i)) {
            glob_reg_state->my_readfrom[i] = MAX(glob_reg_state->my_readfrom[i],
                                                 glob_reg_state->final_setnr);
        }
        if(instr_writes_to_reg(instr, i)) {
            glob_reg_state->my_writtento[i] = MAX(glob_reg_state->my_writtento[i],
                                                  glob_reg_state->final_setnr);
        }
    }
  }
  glob_reg_state->num_sets = MAX(glob_reg_state->num_sets, glob_reg_state->final_setnr);
}

static dr_emit_flags_t
event_basic_block(void *drcontext, void *tag, instrlist_t *bb,
                  bool for_trace, bool translating)
{
    instr_t *instr, *first = instrlist_first(bb);
    uint flags;
    uint cur_flop_count = 0;
#ifdef VERBOSE
    dr_printf("in dynamorio_basic_block(tag="PFX")\n", tag);
# ifdef VERBOSE_VERBOSE
    instrlist_disassemble(drcontext, tag, bb, STDOUT);
# endif
#endif

    /* we use fp ops so we have to save fp state */
    byte fp_raw[512 + 16];
    byte *fp_align = (byte *) ( (((ptr_uint_t)fp_raw) + 16) & ((ptr_uint_t)-16) );


    if (translating) {
        return DR_EMIT_DEFAULT;
    }
    proc_save_fpstate(fp_align);

    int my_readfrom[DR_REG_LAST_VALID_ENUM+MY_NUM_EFLAGS+1];
    int my_writtento[DR_REG_LAST_VALID_ENUM+MY_NUM_EFLAGS+1];
    int i = 0;
    for (i = 0; i < DR_REG_LAST_VALID_ENUM+MY_NUM_EFLAGS+1; i++) {
        my_readfrom[i] = 0;
        my_writtento[i] = 0;
    }

    t_glob_reg_state glob_reg_state = {0,0,0,0,0,0,my_readfrom,my_writtento};

    int my_cur_size = 0;
    for (instr = instrlist_first(bb); instr != NULL; instr = instr_get_next(instr)) {
        my_cur_size++;

        /* ILP Calculations */
        glob_reg_state.raw_setnr = 1;
        glob_reg_state.war_setnr = 1;
        glob_reg_state.waw_setnr = 1;
        glob_reg_state.else_setnr = 1;
        glob_reg_state.final_setnr = 1;
        calc_set_num(instr, &glob_reg_state);

        /* Count flops */
        if( instr_is_floating( instr ) ) {
            cur_flop_count += 1;
        }
    }

    //now we can calculate the ILP.
    float ilp = ((float)my_cur_size) / ((float)(glob_reg_state.num_sets != 0 ?
                glob_reg_state.num_sets : 1));

    dr_mutex_lock(stats_mutex);

    // Due to lack of memory, we only store the ILPs for the latest MY_MAX_BB
    // basic blocks. This enables us to run e.g. firefox.
    int my_cur_num = my_bbcount % MY_MAX_BB;
    my_bbcount++;
    if(my_cur_num == 0 && my_bbcount > 1) {
         dr_printf("Overflow at %d\n", my_bbcount);
    }
    my_bbexecs[my_cur_num] = 0; //initialize
    my_bbsizes[my_cur_num] = my_cur_size;
    bb_flop_count[my_cur_num] = cur_flop_count;
    my_bbilp[my_cur_num] = ilp;

    dr_mutex_unlock(stats_mutex);

#ifdef USE_CLEAN_CALL
     dr_insert_clean_call(drcontext, bb, instrlist_first(bb), clean_call, false, 1,
                           OPND_CREATE_INT32(my_cur_num));
#else
#ifdef INSERT_AT_END
    instr = NULL;
#else
    // Find place to insert inc instruction
    for (instr = first; instr != NULL; instr = instr_get_next(instr)) {
        flags = instr_get_arith_flags(instr);
        if (TESTALL(EFLAGS_WRITE_6, flags) && !TESTANY(EFLAGS_READ_6, flags))
            break;
    }
#endif
    if (instr == NULL) { // no suitable place found, save regs
        dr_save_reg(drcontext, bb, first, DR_REG_XAX, SPILL_SLOT_1);
        dr_save_arith_flags_to_xax(drcontext, bb, first);
    }
    // Increment my_bbexecs[my_current_bb] using the lock prefix
    instrlist_meta_preinsert
        (bb, (instr == NULL) ? first : instr,
         LOCK(INSTR_CREATE_inc(drcontext, OPND_CREATE_ABSMEM
                               ((byte *)&(my_bbexecs[my_cur_num]), OPSZ_4))));
    if (instr == NULL) { // no suitable place found earlier, restore regs
        dr_restore_arith_flags_from_xax(drcontext, bb, first);
        dr_restore_reg(drcontext, bb, first, DR_REG_XAX, SPILL_SLOT_1);
    }
#endif

    proc_restore_fpstate(fp_align);
    
#if defined(VERBOSE) && defined(VERBOSE_VERBOSE)
    dr_printf("Finished instrumenting dynamorio_basic_block(tag="PFX")\n", tag);
    instrlist_disassemble(drcontext, tag, bb, STDOUT);
#endif
    return DR_EMIT_DEFAULT;
}

#ifdef USE_CLEAN_CALL
static void clean_call(uint cur_num){
    dr_mutex_lock(stats_mutex);
    my_bbexecs[cur_num]++;
    dr_mutex_unlock(stats_mutex);
}
#endif

inline void update_setnrs(instr_t *instr, int dr_reg_enum_start, int dr_reg_enum_stop,
                          t_glob_reg_state* glob_reg_state) {
    int i;
    for (i = dr_reg_enum_start; i <= dr_reg_enum_stop; i++) {
        if(i != DR_REG_NULL && i != DR_REG_INVALID) {

            if(instr_reads_from_reg(instr, i)){
                //determine set number for rule 1: RaW (WRITTEN TO + 1)
                glob_reg_state->raw_setnr = MAX(glob_reg_state->raw_setnr,
                                                glob_reg_state->my_writtento[i]+1);
            }

            if(instr_writes_to_reg(instr, i)){
                //determine set number for rule 2: WaR (READ FROM + 1)
                glob_reg_state->war_setnr = MAX(glob_reg_state->war_setnr,
                                                glob_reg_state->my_readfrom[i]+1);

                //determine set number for rule 3: WaW (WRITTEN TO + 1)
                glob_reg_state->waw_setnr = MAX(glob_reg_state->waw_setnr,
                                                glob_reg_state->my_writtento[i]+1);
            }
        }
    }
}

inline void update_eflag_setnrs(instr_t *instr, t_glob_reg_state *glob_reg_state){
    uint flags = instr_get_eflags(instr);
    uint read_masks[MY_NUM_EFLAGS] = {
        EFLAGS_READ_CF,
        EFLAGS_READ_PF,
        EFLAGS_READ_AF,
        EFLAGS_READ_ZF,
        EFLAGS_READ_SF,
        EFLAGS_READ_TF,
        EFLAGS_READ_IF,
        EFLAGS_READ_DF,
        EFLAGS_READ_OF,
        EFLAGS_READ_NT,
        EFLAGS_READ_RF
    };
    uint write_masks[MY_NUM_EFLAGS] = {
        EFLAGS_WRITE_CF,
        EFLAGS_WRITE_PF,
        EFLAGS_WRITE_AF,
        EFLAGS_WRITE_ZF,
        EFLAGS_WRITE_SF,
        EFLAGS_WRITE_TF,
        EFLAGS_WRITE_IF,
        EFLAGS_WRITE_DF,
        EFLAGS_WRITE_OF,
        EFLAGS_WRITE_NT,
        EFLAGS_WRITE_RF
    };

    int i;
    for(i = 0; i < MY_NUM_EFLAGS; i++){
        if(TESTALL(read_masks[i], flags)){
            //determine set number for rule 1: RaW (WRITTEN TO + 1)
            glob_reg_state->raw_setnr = MAX(glob_reg_state->raw_setnr,
                                            glob_reg_state->my_writtento[MY_EFLAGS_OFFSET+i]+1);
        }

        if(TESTALL(write_masks[i], flags)){
            //determine set number for rule 2: WaR (READ FROM + 1)
            glob_reg_state->war_setnr = MAX(glob_reg_state->war_setnr,
                                            glob_reg_state->my_readfrom[MY_EFLAGS_OFFSET+i]+1);

            //determine set number for rule 3: WaW (WRITTEN TO + 1)
            glob_reg_state->waw_setnr = MAX(glob_reg_state->waw_setnr,
                                            glob_reg_state->my_writtento[MY_EFLAGS_OFFSET+i]+1);
        }
    }
}
