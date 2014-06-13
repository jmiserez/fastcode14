#ifndef PYRAMIDS_H
#define PYRAMIDS_H

#include "general_headers.h"

void pyramids(uint32_t nimages, uint32_t nlev, uint32_t r, uint32_t c, double **I, double **W,
              double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize,
              double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
              double ***pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c);
void blend(uint32_t nimages, uint32_t nlev,
           double **pyr, uint32_t *pyr_r, uint32_t *pyr_c,
           double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
           double ***pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c);
void reconstruct_laplacian_pyramid(uint32_t nlev, double *tmp_fullsize, double *tmp2_fullsize, double *tmp_halfsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst);

void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void downsample(double *im, uint32_t r, uint32_t c, double *tmp_halfsize, uint32_t down_r, uint32_t down_c, double *dst);
void downsample_1channel(double *im, uint32_t r, uint32_t c, double *tmp_halfsize, uint32_t down_r, uint32_t down_c, double *dst);
void upsample(double *im, uint32_t r, uint32_t c, uint32_t up_r, uint32_t up_c, double *tmp_halfsize, double *dst);

#endif
