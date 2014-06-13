#include "general_headers.h"
#include "weights_naive_options.h"
#include "basics_naive.h"
#include "pyramids.h"


//
// Exposure Fusion functionality
//

double* fusion_compute(double** I, double contrast_parm, double sat_parm, double wexp_parm,
                        void* _segments){
    segments_t *mem = _segments;

    uint32_t r = mem->r;
    uint32_t c = mem->c;
    uint32_t npixels = mem->npixels;
//    uint32_t nvals = mem->nvals;
    uint32_t nlev = mem->nlev;
    uint32_t nimages = mem->nimages;

    double *R = mem->R;
    double **W = mem->W;
    double **pyr = mem->pyr;
    uint32_t *pyr_r = mem->pyr_r;
    uint32_t *pyr_c = mem->pyr_c;
    double ***pyrW = mem->pyrW;
    uint32_t **pyrW_r = mem->pyrW_r;
    uint32_t **pyrW_c = mem->pyrW_c;
    double ***pyrI = mem->pyrI;
    uint32_t **pyrI_r = mem->pyrI_r;
    uint32_t **pyrI_c = mem->pyrI_c;

    double *tmp_weights = mem->tmp_weights;
    double *tmp2_weights = mem->tmp2_weights;
    double *tmp_fullsize = mem->tmp_fullsize;
    double *tmp2_fullsize = mem->tmp2_fullsize;
    double *tmp_halfsize = mem->tmp_halfsize;
    double *tmp_quartsize = mem->tmp_quartsize;
    double *tmp2_quartsize = mem->tmp2_quartsize;
PERF_FUNC_ENTER
    weights(nimages,npixels,r,c,contrast_parm,sat_parm,wexp_parm,I,W,tmp_weights,tmp2_weights);
PERF_FUNC_EXIT

#ifdef NO_PYRAMIDS
    //for performance measurements we may choose to focus only on weights
    return R;
#endif


    pyramids(nimages, nlev, r, c, I, W,
                  tmp_halfsize, tmp_quartsize, tmp2_quartsize,
                  pyrW, pyrW_r, pyrW_c,
                  pyrI, pyrI_r, pyrI_c);
    blend(nimages, nlev,
               pyr, pyr_r, pyr_c,
               pyrW, pyrW_r, pyrW_c,
               pyrI, pyrI_r, pyrI_c);

    //reconstruct laplacian pyramid
    reconstruct_laplacian_pyramid(nlev,tmp_fullsize,tmp2_fullsize,tmp_halfsize,pyr,pyr_r,pyr_c,r,c,R);
    return R;
}

