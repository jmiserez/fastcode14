#include <stdint.h>

/**
 * Helper macros to compute cost of our computations.
 *
 * Usage:
 *   COST_MEASURE
 *     If defined, counters are activated.
 *     Otherwise, no code is added to source.
 *     NOTE: this must be defined before the header is included
 *
 *   COST_VARIABLES_HERE
 *     Declare this once for your project, e.g. in the main file.
 *
 *   COST_INC_yyy(x)
 *     Increases counter 'yyy' by value 'x'.
 *
 *   COST_yyy
 *     Value of counter 'yyy'.
 *
 *   Valid counters are: ADD, MUL, DIV, POW, ABS, SQRT.
 */

#define COST_T   uint64_t
#define PRI_COST PRId64

#ifndef COST_MEASURE

#define COST_VARIABLES_HERE
#define COST_MODEL_RESET
#define COST_MODEL_LOADF()
#define COST_FUNC_ENTER
#define COST_FUNC_EXIT

#define COST_LOAD ((const COST_T) 0)
#define COST_STORE ((const COST_T) 0)
#define COST_ADD ((const COST_T) 0)
#define COST_MUL ((const COST_T) 0)
#define COST_DIV ((const COST_T) 0)
#define COST_POW ((const COST_T) 0)
#define COST_ABS ((const COST_T) 0)
#define COST_SQRT ((const COST_T) 0)
#define COST_EXP ((const COST_T) 0)

#define COSTF_LOAD ((const COST_T) 0)
#define COSTF_STORE ((const COST_T) 0)
#define COSTF_ADD ((const COST_T) 0)
#define COSTF_MUL ((const COST_T) 0)
#define COSTF_DIV ((const COST_T) 0)
#define COSTF_POW ((const COST_T) 0)
#define COSTF_ABS ((const COST_T) 0)
#define COSTF_SQRT ((const COST_T) 0)
#define COSTF_EXP ((const COST_T) 0)


#define COST_INC_LOAD(x)
#define COST_INC_STORE(x)
#define COST_INC_ADD(x)
#define COST_INC_MUL(x)
#define COST_INC_DIV(x)
#define COST_INC_POW(x)
#define COST_INC_ABS(x)
#define COST_INC_SQRT(x)
#define COST_INC_EXP(x)

#else

extern COST_T __cost_load;
extern COST_T __cost_store;
extern COST_T __cost_add;
extern COST_T __cost_mul;
extern COST_T __cost_div;
extern COST_T __cost_pow;
extern COST_T __cost_abs;
extern COST_T __cost_sqrt;
extern COST_T __cost_exp;

extern COST_T __costF_load;
extern COST_T __costF_store;
extern COST_T __costF_add;
extern COST_T __costF_mul;
extern COST_T __costF_div;
extern COST_T __costF_pow;
extern COST_T __costF_abs;
extern COST_T __costF_sqrt;
extern COST_T __costF_exp;

extern COST_T __costM_load;
extern COST_T __costM_store;
extern COST_T __costM_add;
extern COST_T __costM_mul;
extern COST_T __costM_div;
extern COST_T __costM_pow;
extern COST_T __costM_abs;
extern COST_T __costM_sqrt;
extern COST_T __costM_exp;

#define COST_VARIABLES_HERE  \
    COST_T __costF_load = 0;   \
    COST_T __costF_store = 0;   \
    COST_T __costF_add = 0;   \
    COST_T __costF_mul = 0;   \
    COST_T __costF_div = 0;   \
    COST_T __costF_pow = 0;   \
    COST_T __costF_abs = 0;   \
    COST_T __costF_sqrt = 0;  \
    COST_T __costF_exp = 0;  \
    COST_T __costF_cmp = 0;   \
    COST_T __costM_load = 0;   \
    COST_T __costM_store = 0;   \
    COST_T __costM_add = 0;   \
    COST_T __costM_mul = 0;   \
    COST_T __costM_div = 0;   \
    COST_T __costM_pow = 0;   \
    COST_T __costM_abs = 0;   \
    COST_T __costM_sqrt = 0;  \
    COST_T __costM_exp = 0;  \
    COST_T __costM_cmp = 0;   \
    COST_T __cost_load = 0;   \
    COST_T __cost_store = 0;   \
    COST_T __cost_add = 0;   \
    COST_T __cost_mul = 0;   \
    COST_T __cost_div = 0;   \
    COST_T __cost_pow = 0;   \
    COST_T __cost_abs = 0;   \
    COST_T __cost_sqrt = 0;  \
    COST_T __cost_exp = 0;

#define COST_MODEL_RESET \
    __costF_load = 0;   \
    __costF_store = 0;   \
    __costF_add = 0;   \
    __costF_mul = 0;   \
    __costF_div = 0;   \
    __costF_pow = 0;   \
    __costF_abs = 0;   \
    __costF_sqrt = 0;  \
    __costF_exp = 0;  \
    __costF_cmp = 0;   \
    __costM_load = 0;   \
    __costM_store = 0;   \
    __costM_add = 0;   \
    __costM_mul = 0;   \
    __costM_div = 0;   \
    __costM_pow = 0;   \
    __costM_abs = 0;   \
    __costM_sqrt = 0;  \
    __costM_exp = 0;  \
    __costM_cmp = 0;   \
    __cost_load = 0;      \
    __cost_store = 0;      \
    __cost_add = 0;      \
    __cost_mul = 0;      \
    __cost_div = 0;      \
    __cost_pow = 0;      \
    __cost_abs = 0;      \
    __cost_sqrt = 0;      \
    __cost_exp = 0;

#define COST_MODEL_LOADF { __cost_load = __costF_load; \
    __cost_store = __costF_store; \
    __cost_add = __costF_add; \
    __cost_mul = __costF_mul; \
    __cost_div = __costF_div; \
    __cost_pow = __costF_pow; \
    __cost_abs = __costF_abs; \
    __cost_sqrt = __costF_sqrt; \
    __cost_exp = __costF_exp; }

#define COST_LOAD ((const COST_T) __cost_load)
#define COST_STORE ((const COST_T) __cost_store)
#define COST_ADD ((const COST_T) __cost_add)
#define COST_MUL ((const COST_T) __cost_mul)
#define COST_DIV ((const COST_T) __cost_div)
#define COST_POW ((const COST_T) __cost_pow)
#define COST_ABS ((const COST_T) __cost_abs)
#define COST_SQRT ((const COST_T) __cost_sqrt)
#define COST_EXP ((const COST_T) __cost_exp)

#define COSTF_LOAD ((const COST_T) __costF_load)
#define COSTF_STORE ((const COST_T) __costF_store)
#define COSTF_ADD ((const COST_T) __costF_add)
#define COSTF_MUL ((const COST_T) __costF_mul)
#define COSTF_DIV ((const COST_T) __costF_div)
#define COSTF_POW ((const COST_T) __costF_pow)
#define COSTF_ABS ((const COST_T) __costF_abs)
#define COSTF_SQRT ((const COST_T) __costF_sqrt)
#define COSTF_EXP ((const COST_T) __costF_exp)

#define COST_INC_LOAD(x) {__cost_load += (x);}
#define COST_INC_STORE(x) {__cost_store += (x);}
#define COST_INC_ADD(x) {__cost_add += (x);}
#define COST_INC_MUL(x) {__cost_mul += (x);}
#define COST_INC_DIV(x) {__cost_div += (x);}
#define COST_INC_POW(x) {__cost_pow += (x);}
#define COST_INC_ABS(x) {__cost_abs += (x);}
#define COST_INC_SQRT(x) {__cost_sqrt += (x);}
#define COST_INC_EXP(x) {__cost_exp += (x);}

#define COST_FUNC_ENTER {__costM_load = __cost_load; \
                         __costM_store = __cost_store; \
                         __costM_add = __cost_add; \
                         __costM_mul = __cost_mul; \
                         __costM_div = __cost_div; \
                         __costM_pow = __cost_pow; \
                         __costM_abs = __cost_abs; \
                         __costM_sqrt = __cost_sqrt; \
                         __costM_exp = __cost_exp;}

#define COST_FUNC_EXIT {__costF_load += __cost_load-__costM_load; \
                         __costF_store += __cost_store - __costM_store; \
                         __costF_add += __cost_add - __costM_add; \
                         __costF_mul += __cost_mul - __costM_mul; \
                         __costF_div += __cost_div - __costM_div; \
                         __costF_pow += __cost_pow - __costM_pow; \
                         __costF_abs += __cost_abs - __costM_abs; \
                         __costF_sqrt += __cost_sqrt - __costM_sqrt; \
                         __costF_exp += __cost_exp - __costM_exp;}

#endif
