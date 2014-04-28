#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include<stddef.h>

void conv3x3_monochrome_replicate( double* dst, double* im, int w, int h,
                                   double* f );

#endif // CONVOLUTION_H
