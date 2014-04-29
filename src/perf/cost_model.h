#include <stdint.h>

#define COST_T   uint64_t
#define PRI_COST PRId64

#ifndef COST_MEASURE

#define COST_VARIABLES_HERE

#define COST_ADD ((const COST_T) 0)
#define COST_MUL ((const COST_T) 0)
#define COST_DIV ((const COST_T) 0)
#define COST_POW ((const COST_T) 0)
#define COST_ABS ((const COST_T) 0)

#define COST_INC_ADD(x)
#define COST_INC_MUL(x)
#define COST_INC_DIV(x)
#define COST_INC_POW(x)
#define COST_INC_ABS(x)

#else

extern COST_T __cost_add;
extern COST_T __cost_mul;
extern COST_T __cost_div;
extern COST_T __cost_pow;
extern COST_T __cost_abs;

#define COST_VARIABLES_HERE  \
    COST_T __cost_add = 0;   \
    COST_T __cost_mul = 0;   \
    COST_T __cost_div = 0;   \
    COST_T __cost_pow = 0;   \
    COST_T __cost_abs = 0;


#define COST_ADD ((const COST_T) __cost_add)
#define COST_MUL ((const COST_T) __cost_mul)
#define COST_DIV ((const COST_T) __cost_div)
#define COST_POW ((const COST_T) __cost_pow)
#define COST_ABS ((const COST_T) __cost_abs)

#define COST_INC_ADD(x) {__cost_add += (x);}
#define COST_INC_MUL(x) {__cost_mul += (x);}
#define COST_INC_DIV(x) {__cost_div += (x);}
#define COST_INC_POW(x) {__cost_pow += (x);}
#define COST_INC_ABS(x) {__cost_abs += (x);}

#endif
