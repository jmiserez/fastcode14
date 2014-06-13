#ifndef GENERAL_HEADERS_H
#define GENERAL_HEADERS_H

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

// When using this, the -Winline compiler flag needs to be set.
// Then an error will occur when a function was not inlined.
#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE
#endif

//104x104 pixels fit in L2 cache (256KB)
//512x512 pixels fit in L3 cache (6MB)
#define L1_CACHE_SIZE ((L1_CACHE_KB * 1024) / sizeof(double))
#define L2_CACHE_SIZE ((L2_CACHE_KB * 1024) / sizeof(double))
#define LLC_CACHE_SIZE ((LLC_CACHE_KB * 1024) / sizeof(double))
#define CACHE_LINE_SIZE (CACHE_LINE_BYTES / sizeof(double))

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define CHANNELS ( 3 )

//NOTE: image size MUST be an exact multiple of blocksize
#define NB 128

#endif
