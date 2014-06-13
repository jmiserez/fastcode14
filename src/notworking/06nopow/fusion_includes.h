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

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE
#endif

#define L1_CACHE_SIZE ((L1_CACHE_KB * 1024) / sizeof(double))
#define L2_CACHE_SIZE ((L2_CACHE_KB * 1024) / sizeof(double))
#define LLC_CACHE_SIZE ((LLC_CACHE_KB * 1024) / sizeof(double))
#define CACHE_LINE_SIZE (CACHE_LINE_BYTES / sizeof(double))

#define TMP_WEIGHTS_SIZE 512*512*3

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define CHANNELS 3

#endif // FUSION_INCLUDES_H
