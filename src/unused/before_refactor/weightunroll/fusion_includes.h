#ifndef FUSION_INCLUDES_H
#define FUSION_INCLUDES_H

#include <stdio.h>
#include <string.h>
#include <stdint.h>
//#include <stdint-gcc.h> Not available on RedHat6@CAB. Not needed.
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <tiffio.h>
#include "fusion.h"
#include "fusion_perfcost.h"
#include <x86intrin.h>

#define FORCE_INLINE __attribute__((always_inline)) inline

#define L1_CACHE_SIZE ((L1_CACHE_KB * 1024) / sizeof(double))
#define L2_CACHE_SIZE ((L2_CACHE_KB * 1024) / sizeof(double))
#define LLC_CACHE_SIZE ((LLC_CACHE_KB * 1024) / sizeof(double))
#define CACHE_LINE_SIZE (CACHE_LINE_BYTES / sizeof(double))

#define CHANNELS ( 3 )


#define WEIGHTS_UNROLL_NB_R ( 1 )
#define WEIGHTS_UNROLL_NB_C ( 2 )


#define WEIGHTS_ELEMS_NB ( L2_CACHE_SIZE/CHANNELS )
#define WEIGHTS_SQRT_NB ( (int)floor(sqrt( WEIGHTS_ELEMS_NB )) )
#define WEIGHTS_MIN_NB ( WEIGHTS_SQRT_NB - (WEIGHTS_SQRT_NB % CACHE_LINE_SIZE) )
#define WEIGHTS_MAX_NB ( WEIGHTS_ELEMS_NB/WEIGHTS_MIN_NB )

//#define WEIGHTS_NB_R ( WEIGHTS_MAX_NB - WEIGHTS_UNROLL_NB_R - (WEIGHTS_MAX_NB % WEIGHTS_UNROLL_NB_R) + 2)
//#define WEIGHTS_NB_C ( WEIGHTS_MIN_NB - WEIGHTS_UNROLL_NB_C - (WEIGHTS_MIN_NB % WEIGHTS_UNROLL_NB_C) + 2)

#define WEIGHTS_NB_R ( 32 )
#define WEIGHTS_NB_C ( 32 )


#define TMP_WEIGHTS_SIZE ( (WEIGHTS_NB_R+2)*(WEIGHTS_NB_C+2)*CHANNELS )


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#endif // FUSION_INCLUDES_H
