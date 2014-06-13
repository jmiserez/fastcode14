#ifndef WEIGHTS_H
#define WEIGHTS_H

#include "general_headers.h"

void weights(uint32_t nimages, uint32_t npixels, uint32_t r, uint32_t c,
             double **I, double **W, double *tmp_weights, double *tmp2_weights);
void conv3x3_mult_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* factors, double* dst);
#endif
