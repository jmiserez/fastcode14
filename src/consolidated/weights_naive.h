#ifndef WEIGHTS_H
#define WEIGHTS_H

#include "general_headers.h"

void weights(uint32_t nimages, uint32_t npixels, uint32_t r, uint32_t c,
             double **I, double **W, double *tmp_weights, double *tmp2_weights);
void rgb2gray(double *im, size_t npixels, double* dst);
void conv3x3_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* dst);

#endif
