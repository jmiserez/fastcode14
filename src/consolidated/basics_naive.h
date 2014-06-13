#ifndef BASICS_H
#define BASICS_H

#include "general_headers.h"

int malloc_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr, uint32_t **pyr_r, uint32_t **pyr_c, bool level0_is_ref);
void free_pyramid(uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, bool level0_is_ref);

typedef struct {
    uint32_t r;
    uint32_t c;
    uint32_t npixels;
    uint32_t nvals;
    uint32_t nlev;
    uint32_t nimages;
    uint32_t max_upsampled_r;
    uint32_t max_upsampled_c;

    double *R;
    double **W;

    //result pyramid
    double **pyr;
    uint32_t *pyr_r;
    uint32_t *pyr_c;

    //pyrW[n] is pyramid of weight map for each of the N images
    double ***pyrW;
    uint32_t **pyrW_r;
    uint32_t **pyrW_c;

    //pyrI[n] is image pyramid for each of the N images
    double ***pyrI;
    uint32_t **pyrI_r;
    uint32_t **pyrI_c;

    //regular size scratch spaces
    double *tmp_weights;
    double *tmp2_weights;
    double *tmp_fullsize;
    double *tmp2_fullsize;
    double *tmp_halfsize;
    double *tmp_quartsize;
    double *tmp2_quartsize;
} segments_t;

#endif
