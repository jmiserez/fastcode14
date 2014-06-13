// matlab-based implementation as shown in meeting
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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define CHANNELS 3

typedef struct {
    uint32_t r;
    uint32_t c;
    uint32_t npixels;
    uint32_t nvals;
    uint32_t nlev;
    uint32_t nimages;
    uint32_t max_upsampled_r;
    uint32_t max_upsampled_c;

    double ***C;
    double *R;
    double **W;

    //result pyramids per channel
    double ***pyr;
    uint32_t *pyr_r;
    uint32_t *pyr_c;

    //pyrW[n] is pyramid of weight map for each of the N images
    double ***pyrW;
    uint32_t **pyrW_r;
    uint32_t **pyrW_c;

    //pyrI[channel][n] is image pyramid for each of the N images, per channel
    double ****pyrI;
    uint32_t **pyrI_r;
    uint32_t **pyrI_c;

    //regular size scratch spaces
    double *tmp_fullsize;
    double *tmp2_fullsize;
    double *tmp_halfsize;
    double *tmp_quartsize;
    double *tmp2_quartsize;
} segments_t;


int malloc_pyramid(uint32_t r, uint32_t c, uint32_t nlev, double ***pyr, uint32_t *pyr_r, uint32_t *pyr_c, bool level0_is_ref);
void free_pyramid(uint32_t nlev, double **pyr, bool level0_is_ref);

void weights(uint32_t nimages, uint32_t npixels, uint32_t r, uint32_t c, double contrast_parm, double sat_parm, double wexp_parm, double **I, double ***C, double **W);
void pyramids(uint32_t nimages, uint32_t nlev, uint32_t r, uint32_t c, double ***C, double **W,
              double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize,
              double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
              double ****pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c);
void blend(uint32_t nimages, uint32_t nlev,
           double ***pyr, uint32_t *pyr_r, uint32_t *pyr_c,
           double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
           double ****pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c);
void reconstruct_laplacian_pyramid(uint32_t nlev, uint32_t npixels, double *tmp_fullsize, double *tmp2_fullsize, double *tmp_halfsize, double ***pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst);
void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void downsample(double *im, uint32_t r, uint32_t c, double *tmp_halfsize, uint32_t down_r, uint32_t down_c, double *dst);
void upsample(double *im, uint32_t r, uint32_t c, uint32_t up_r, uint32_t up_c, double *tmp_halfsize, double *dst);



//
// Memory handling
//

int fusion_alloc(void** _segments, int w, int h, int N){

    segments_t *mem = malloc(sizeof(segments_t));
    if(mem == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->npixels = w*h;
    mem->nvals = w*h*3;
    mem->nlev = (uint32_t)(floor((log2(MIN(w,h)))));
    mem->nimages = N;
    mem->r = h;
    mem->c = w;
    mem->max_upsampled_r = h + (h%2);
    mem->max_upsampled_c = w + (w%2);

    mem->C = calloc(mem->nimages,sizeof(double**));
    if(mem->C == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    for (int n = 0; n < mem->nimages; n++){
        mem->C[n] = calloc(CHANNELS,sizeof(double*));
        if(mem->C[n] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
        for (int k = 0; k < CHANNELS; k++){
            mem->C[n][k] = calloc(mem->npixels,sizeof(double));
            if(mem->C[n][k] == NULL){
                return FUSION_ALLOC_FAILURE;
            }
        }
    }

    mem->R = calloc(mem->nvals,sizeof(double));
    if(mem->R == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->W = calloc(mem->nimages,sizeof(double*));
    if(mem->W == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    for (int n = 0; n < mem->nimages; n++){
        mem->W[n] = calloc(mem->npixels,sizeof(double));
        if(mem->W[n] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
    }

    mem->pyr = calloc(CHANNELS,sizeof(double**));
    if(mem->pyr == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->pyrW = calloc(mem->nimages,sizeof(double**));
    if(mem->pyrW == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrW_r = calloc(mem->nimages,sizeof(double*));
    if(mem->pyrW_r == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrW_c = calloc(mem->nimages,sizeof(double*));
    if(mem->pyrW_c == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->pyrI = calloc(CHANNELS,sizeof(double***));
    if(mem->pyrI == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    for(int k = 0; k < CHANNELS; k++){
        mem->pyrI[k] = calloc(mem->nimages,sizeof(double**));
        if(mem->pyrI[k] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
    }
    mem->pyrI_r = calloc(mem->nimages,sizeof(double*));
    if(mem->pyrI_r == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrI_c = calloc(mem->nimages,sizeof(double*));
    if(mem->pyrI_c == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    //alloc pyramids

    for (int n = 0; n < N; n++){
        mem->pyrW_r[n] = (uint32_t*) calloc(mem->nlev,sizeof(uint32_t));
        if(mem->pyrW_r[n] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
        mem->pyrW_c[n] = (uint32_t*) calloc(mem->nlev,sizeof(uint32_t));
        if(mem->pyrW_c[n] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
        if(malloc_pyramid(mem->r,mem->c,mem->nlev,&(mem->pyrW[n]), mem->pyrW_r[n], mem->pyrW_c[n],true) != FUSION_ALLOC_SUCCESS){
            return FUSION_ALLOC_FAILURE;
        }
    }

    mem->pyr_r = (uint32_t*) calloc(mem->nlev,sizeof(uint32_t));
    if(mem->pyr_r == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyr_c = (uint32_t*) calloc(mem->nlev,sizeof(uint32_t));
    if(mem->pyr_c == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    for(int k = 0; k < CHANNELS; k++){
        if(malloc_pyramid(mem->r,mem->c,mem->nlev,&(mem->pyr[k]), mem->pyr_r, mem->pyr_c,false) != FUSION_ALLOC_SUCCESS){
            return FUSION_ALLOC_FAILURE;
        }
    }

    for (int n = 0; n < N; n++){
        mem->pyrI_r[n] = (uint32_t*) calloc(mem->nlev,sizeof(uint32_t));
        if(mem->pyrI_r[n] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
        mem->pyrI_c[n] = (uint32_t*) calloc(mem->nlev,sizeof(uint32_t));
        if(mem->pyrI_c[n] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
        for(int k = 0; k < CHANNELS; k++){
            if(malloc_pyramid(mem->r,mem->c,mem->nlev,&(mem->pyrI[k][n]), mem->pyrI_r[n], mem->pyrI_c[n],false) != FUSION_ALLOC_SUCCESS){
                return FUSION_ALLOC_FAILURE;
            }
        }
    }

    //alloc tmp vars

    mem->tmp_fullsize = calloc(mem->max_upsampled_r * mem->max_upsampled_c,sizeof(double));
    if(mem->tmp_fullsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp2_fullsize = calloc(mem->max_upsampled_r * mem->max_upsampled_c,sizeof(double));
    if(mem->tmp2_fullsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp_halfsize = calloc(
                (mem->max_upsampled_r * mem->max_upsampled_c)+1 / 2, sizeof(double));
    if(mem->tmp_halfsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp_quartsize = calloc(
                ((mem->max_upsampled_r * mem->max_upsampled_c)+1 / 2)+1 / 2, sizeof(double));
    if(mem->tmp_quartsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp2_quartsize = calloc(
                ((mem->max_upsampled_r * mem->max_upsampled_c)+1 / 2)+1 / 2, sizeof(double));
    if(mem->tmp2_quartsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    *_segments = mem;

    return FUSION_ALLOC_SUCCESS;
}

/**
 * Allocate memory for the gaussian/laplacian pyramid at *pyr
 */
int malloc_pyramid(uint32_t r, uint32_t c, uint32_t nlev, double ***pyr, uint32_t *pyr_r, uint32_t *pyr_c, bool level0_is_ref){
    *pyr = (double**) calloc(nlev,sizeof(double*));
    assert(*pyr != NULL);

    uint32_t r_level = r;
    uint32_t c_level = c;

    for(int i = 0; i < nlev; i++){
        pyr_r[i] = r_level; //store dimension r at level i
        pyr_c[i] = c_level; //store dimension c at level i

        if(i != 0 || !level0_is_ref){
            size_t L_len = r_level*c_level;
            double* L = calloc(L_len,sizeof(double));
            if(L == NULL){
                return FUSION_ALLOC_FAILURE;
            }
            (*pyr)[i] = L; //add entry to array of pointers to image levels
        } else {
            (*pyr)[i] = NULL;
        }
        // for next level
        r_level = r_level / 2 + (r_level % 2);
        c_level = c_level / 2 + (c_level % 2);
    }
    return FUSION_ALLOC_SUCCESS;
}

void fusion_free( void* _segments ){
    segments_t *mem = _segments;

    free(mem->tmp_fullsize);
    free(mem->tmp2_fullsize);
    free(mem->tmp_halfsize);
    free(mem->tmp_quartsize);
    free(mem->tmp2_quartsize);

    for(int k = 0; k < CHANNELS; k++){
        free_pyramid(mem->nlev,mem->pyr[k],false);
    }
    free(mem->pyr_r);
    free(mem->pyr_c);

    for (int n = 0; n < mem->nimages; n++){
        free(mem->W[n]);
        free_pyramid(mem->nlev,mem->pyrW[n], true);
        for(int k = 0; k < CHANNELS; k++){
            free_pyramid(mem->nlev,mem->pyrI[k][n], false);
        }
        free(mem->pyrI_r[n]);
        free(mem->pyrI_c[n]);
    }


    free(mem->pyrI);
    free(mem->pyrI_r);
    free(mem->pyrI_c);
    free(mem->pyrW);
    free(mem->pyrW_r);
    free(mem->pyrW_c);
    free(mem->W);

    for (int n = 0; n < mem->nimages; n++){
        for(int k = 0; k < CHANNELS; k++){
            free(mem->C[n][k]);
        }
        free(mem->C[n]);
    }
    free(mem->C);
    free(mem->R);
    free(_segments);
}

void free_pyramid(uint32_t nlev, double **pyr, bool level0_is_ref){
    for(int i = 0; i < nlev; i++){
        if(i != 0 || !level0_is_ref){
            free(pyr[i]);
        }
    }
    free(pyr);
}

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
    double ***C = mem->C;
    double **W = mem->W;
    double ***pyr = mem->pyr;
    uint32_t *pyr_r = mem->pyr_r;
    uint32_t *pyr_c = mem->pyr_c;
    double ***pyrW = mem->pyrW;
    uint32_t **pyrW_r = mem->pyrW_r;
    uint32_t **pyrW_c = mem->pyrW_c;
    double ****pyrI = mem->pyrI;
    uint32_t **pyrI_r = mem->pyrI_r;
    uint32_t **pyrI_c = mem->pyrI_c;

    double *tmp_fullsize = mem->tmp_fullsize;
    double *tmp2_fullsize = mem->tmp2_fullsize;
    double *tmp_halfsize = mem->tmp_halfsize;
    double *tmp_quartsize = mem->tmp_quartsize;
    double *tmp2_quartsize = mem->tmp2_quartsize;

    PERF_FUNC_ENTER
    weights(nimages,npixels,r,c,contrast_parm,sat_parm,wexp_parm,I,C,W);
    PERF_FUNC_EXIT

    pyramids(nimages, nlev, r, c, C, W,
                  tmp_halfsize, tmp_quartsize, tmp2_quartsize,
                  pyrW, pyrW_r, pyrW_c,
                  pyrI, pyrI_r, pyrI_c);
    blend(nimages, nlev,
               pyr, pyr_r, pyr_c,
               pyrW, pyrW_r, pyrW_c,
               pyrI, pyrI_r, pyrI_c);


    reconstruct_laplacian_pyramid(nlev,npixels,tmp_fullsize,tmp2_fullsize,tmp_halfsize,pyr,pyr_r,pyr_c,r,c,R);
    return R;
}

FORCE_INLINE void _weights_inner(double *im, double **C, int at, int top, int left, int right, int bottom, double *Wn, double sat_parm, double contrast_parm, double wexp_parm){
    int at3 = at*CHANNELS;
    double r = im[at3];
    double g = im[at3+1];
    double b = im[at3+2];

    //store individual channels
    //TODO SIMD nocache
    C[0][at] = r;
    C[1][at] = g;
    C[2][at] = b;

    double r1 = im[top];
    double g1 = im[top+1];
    double b1 = im[top+2];
    double r2 = im[left];
    double g2 = im[left+1];
    double b2 = im[left+2];
    double r3 = im[right];
    double g3 = im[right+1];
    double b3 = im[right+2];
    double r4 = im[bottom];
    double g4 = im[bottom+1];
    double b4 = im[bottom+2];
    double rw = 0.2989;
    double gw = 0.5870;
    double bw = 0.1140;

    double grey = rw * r + gw * g + bw * b;
    double grey1 = rw * r1 + gw * g1 + bw * b1;
    double grey2 = rw * r2 + gw * g2 + bw * b2;
    double grey3 = rw * r3 + gw * g3 + bw * b3;
    double grey4 = rw * r4 + gw * g4 + bw * b4;
    COST_INC_ADD(2);
    COST_INC_MUL(3);
    double mu = (r + g + b) / 3.0;
    COST_INC_ADD(2);
    COST_INC_DIV(1);
    double rmu = r-mu;
    double gmu = g-mu;
    double bmu = b-mu;
    COST_INC_ADD(3);
    double rz = r-0.5;
    double gz = g-0.5;
    double bz = b-0.5;
    COST_INC_ADD(3);
    double rzrz = rz*rz;
    double gzgz = gz*gz;
    double bzbz = bz*bz;
    COST_INC_MUL(3);
    double re = exp(-12.5*rzrz);
    double ge = exp(-12.5*gzgz);
    double be = exp(-12.5*bzbz);
    COST_INC_EXP(3);
    COST_INC_MUL(3);
    double t1 = sqrt((rmu*rmu + gmu*gmu + bmu*bmu)/3.0);
    COST_INC_SQRT(1);
    COST_INC_ADD(3);
    COST_INC_MUL(3);
    COST_INC_DIV(1);
    double t2 = pow(fabs(t1),sat_parm);
    COST_INC_POW(1);
    COST_INC_ABS(1);
    double t3 = re*ge*be;
    COST_INC_MUL(2);
    double t4 = pow(fabs(t3),wexp_parm);
    COST_INC_POW(1);
    COST_INC_ABS(1);

    double t5 = t2 * t4;
    COST_INC_MUL(1);

    double t6 = -4.0*grey+grey1+grey2+grey3+grey4;
    COST_INC_MUL(1);
    COST_INC_ADD(4);

    double t7 = pow(fabs(t6),contrast_parm);
    COST_INC_POW(1);
    COST_INC_ABS(1);

    double t8 = t5 * t7;
    COST_INC_MUL(1);

    double t9 = t8 + 1.0E-12;
    COST_INC_ADD(1);

    Wn[at] = t9;
}

FORCE_INLINE void _weights_inner_row(double *im, double **C, int c, int i, int itop, int ibottom, double *Wn, double sat_parm, double contrast_parm, double wexp_parm){

    //center cols
    for(int j = 1; j < c-1; j++){
        int at = (i*c+j);
        int top = ((itop)*c+j)*CHANNELS;
        int bottom = ((ibottom)*c+j)*CHANNELS;
        int left = (i*c+j-1)*CHANNELS;
        int right = (i*c+j+1)*CHANNELS;
        _weights_inner(im,C,at,top,left,right,bottom,Wn,sat_parm,contrast_parm,wexp_parm);
    }
    //left cols
    int j = 0;
    int at = (i*c+j);
    int top = ((itop)*c+j)*CHANNELS;
    int bottom = ((ibottom)*c+j)*CHANNELS;
    int left = (i*c+j)*CHANNELS;
    int right = (i*c+j+1)*CHANNELS;
    _weights_inner(im,C,at,top,left,right,bottom,Wn,sat_parm,contrast_parm,wexp_parm);
    //right cols
    j = c-1;
    at = (i*c+j);
    top = ((itop)*c+j)*CHANNELS;
    bottom = ((ibottom)*c+j)*CHANNELS;
    left = (i*c+j-1)*CHANNELS;
    right = (i*c+j)*CHANNELS;
    _weights_inner(im,C,at,top,left,right,bottom,Wn,sat_parm,contrast_parm,wexp_parm);
}

void weights(uint32_t nimages, uint32_t npixels, uint32_t r, uint32_t c, double contrast_parm, double sat_parm, double wexp_parm,
             double **I, double ***C, double **W){

    for (int n = 0; n < nimages; n++){
        double *im = I[n];
        double **imC = C[n];
        double *Wn = W[n];

        //center rows
        for(int i = 1; i < r-1; i++){
            int itop = i-1;
            int ibottom = i+1;
            _weights_inner_row(im,imC,c,i,itop,ibottom,Wn,sat_parm,contrast_parm,wexp_parm);
        }
        //top row
        int i = 0;
        int itop = 0;
        int ibottom = 1;
        _weights_inner_row(im,imC,c,i,itop,ibottom,Wn,sat_parm,contrast_parm,wexp_parm);
        //bottom row
        i = r-1;
        itop = r-2;
        ibottom = r-1;
        _weights_inner_row(im,imC,c,i,itop,ibottom,Wn,sat_parm,contrast_parm,wexp_parm);
    }

    //normalization
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            double sum = 0.0; //sum of all weights for this pixel
            for (int n = 0; n < nimages; n++){
                int at = (i*c+j);
                double *Wn = W[n];
                sum += Wn[at]; COST_INC_ADD(1);
            }
            for (int n = 0; n < nimages; n++){
                int at = (i*c+j);
                double *Wn = W[n];
                Wn[at] = Wn[at] / sum; COST_INC_DIV(1); //beware of division by zero
            }
        }
    }
}

void pyramids(uint32_t nimages, uint32_t nlev, uint32_t r, uint32_t c, double ***C, double **W,
              double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize,
              double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
              double ****pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c){
    //multiresolution blending
    for (int n = 0; n < nimages; n++){
        //construct 1-channel gaussian pyramid from weights
        gaussian_pyramid(W[n],r,c,nlev,tmp_halfsize,pyrW[n],pyrW_r[n],pyrW_c[n]);
    }
    for (int n = 0; n < nimages; n++){
        for(int k = 0; k < CHANNELS; k++){
            double ***pyrIk = pyrI[k];
            laplacian_pyramid(C[n][k],r,c,nlev,tmp_halfsize,tmp_quartsize,tmp2_quartsize,pyrIk[n],pyrI_r[n],pyrI_c[n]);
        }
    }
}

void blend(uint32_t nimages, uint32_t nlev,
           double ***pyr, uint32_t *pyr_r, uint32_t *pyr_c,
           double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
           double ****pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c){
    //weighted blend
    for(int k = 0; k < CHANNELS; k++){
        double **pyrk = pyr[k];
        double ***pyrIk = pyrI[k];
        if(0 < nimages){
            int n = 0;
            for(int v = 0; v < nlev; v++){
                for(int i = 0; i < pyrI_r[n][v]; i++){
                    for(int j = 0; j < pyrI_c[n][v]; j++){
                         pyrk[v][(i*pyr_c[v]+j)] = pyrW[n][v][i*pyrI_c[n][v]+j] * pyrIk[n][v][(i*pyrI_c[n][v]+j)];
                    }
                }
            }
        }
        for (int n = 1; n < nimages; n++){
            for(int v = 0; v < nlev; v++){
                for(int i = 0; i < pyrI_r[n][v]; i++){
                    for(int j = 0; j < pyrI_c[n][v]; j++){
                        pyrk[v][(i*pyr_c[v]+j)] += pyrW[n][v][i*pyrI_c[n][v]+j] * pyrIk[n][v][(i*pyrI_c[n][v]+j)];
                    }
                }
            }
        }
    }
}

void reconstruct_laplacian_pyramid(uint32_t nlev, uint32_t npixels, double *tmp_fullsize, double *tmp2_fullsize, double *tmp_halfsize, double ***pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst){

    for(int k = 0; k < CHANNELS; k++){
        double **pyrk = pyr[k];

        // reconstruct image for channel k
        if (nlev-2 >= 0){
            int v = nlev-2;
            upsample(pyrk[v+1],pyr_r[v+1],pyr_c[v+1],pyr_r[v],pyr_c[v],tmp_halfsize,tmp2_fullsize);
            for(int i = 0; i < pyr_r[v]*pyr_c[v]; i++){
                tmp_fullsize[i] = pyrk[v][i] + tmp2_fullsize[i]; COST_INC_ADD(1);
            }
        }
        for (int v = nlev-3; v >= 0; v--){
            upsample(tmp_fullsize,pyr_r[v+1],pyr_c[v+1],pyr_r[v],pyr_c[v],tmp_halfsize,tmp2_fullsize);
            for(int i = 0; i < pyr_r[v]*pyr_c[v]; i++){
                tmp_fullsize[i] = pyrk[v][i] + tmp2_fullsize[i]; COST_INC_ADD(1);
            }
        }

        //store into result image
        //TODO SIMD nocache
        for(int i = 0; i < npixels; i++){
            dst[i*CHANNELS+k] = tmp_fullsize[i];
        }
    }
}

void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    pyr[0] = im;
    if(1 < nlev){
        int v = 1;
        downsample(pyr[0],pyr_r[v-1],pyr_c[v-1],tmp_halfsize,pyr_r[v],pyr_c[v],pyr[v]);
    }
    for(int v = 2; v < nlev; v++){
        //downsample image and store into level
        downsample(pyr[v-1],pyr_r[v-1],pyr_c[v-1],tmp_halfsize,pyr_r[v],pyr_c[v],pyr[v]);
    }
}

void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    uint32_t S_r = r;
    uint32_t S_c = c;
    uint32_t T_r = r;
    uint32_t T_c = c;

    double *tmp = NULL;

    if(0 < nlev-1){
        int v = 0;
        S_r = pyr_r[v+1];
        S_c = pyr_c[v+1];
        downsample(im,T_r,T_c,tmp_halfsize,S_r,S_c,tmp2_quartsize);
        upsample(tmp2_quartsize,S_r,S_c,pyr_r[v],pyr_c[v],tmp_halfsize,pyr[v]);
        for(int i = 0; i < T_r*T_c; i++){
            pyr[v][i] = im[i] - pyr[v][i]; COST_INC_ADD(1);
        }
        T_r = S_r;
        T_c = S_c;
        double *tmp = tmp_quartsize;
        tmp_quartsize = tmp2_quartsize;
        tmp2_quartsize = tmp;
    }
    for(int v = 1; v < nlev-1; v++){
        S_r = pyr_r[v+1];
        S_c = pyr_c[v+1];
        downsample(tmp_quartsize,T_r,T_c,tmp_halfsize,S_r,S_c,tmp2_quartsize);
        upsample(tmp2_quartsize,S_r,S_c,pyr_r[v],pyr_c[v],tmp_halfsize,pyr[v]);
        for(int i = 0; i < T_r*T_c; i++){
            pyr[v][i] = tmp_quartsize[i] - pyr[v][i]; COST_INC_ADD(1);
        }
        T_r = S_r;
        T_c = S_c;
        tmp = tmp_quartsize;
        tmp_quartsize = tmp2_quartsize;
        tmp2_quartsize = tmp;
    }
    //memcpy(dst,src,src_len*sizeof(double));
    for(int i = 0; i < T_r*T_c; i++){
        pyr[nlev-1][i] = tmp_quartsize[i];
    }
}

void downsample(double *im, uint32_t r, uint32_t c, double *tmp_halfsize, uint32_t down_r, uint32_t down_c, double *dst){
    int c2 = c/2+((c-2) % 2); //tmp_halfsize is only half the size
    for(int i = 0; i < r; i++){
        int j_half = 1;
        for(int j = 2; j < c-2; j+=2){ //every 2nd column
            tmp_halfsize[(i*c2+j_half)] =
                    im[((i  )*c+(j-2))]*.0625 +
                    im[((i  )*c+(j-1))]*.25 +
                    im[((i  )*c+(j  ))]*.375 +
                    im[((i  )*c+(j+1))]*.25 +
                    im[((i  )*c+(j+2))]*.0625;
            COST_INC_ADD(4);
            COST_INC_MUL(5);
            j_half++;
        }
        //left edge
        int j = 0; // 1 0 [0 1 2 ... ]
        tmp_halfsize[(i*c2)] =
                im[((i  )*c+(j+1))]*.3125 +
                im[((i  )*c+(j  ))]*.625 +
                im[((i  )*c+(j+2))]*.0625;
        COST_INC_ADD(2);
        COST_INC_MUL(3);
        //right edge
        if((c-2) % 2 == 0){
            j = c-2; // [ ... -2 -1 0 1] 1
            tmp_halfsize[(i*c2+(j/2))] =
                    im[((i  )*c+(j-2))]*.0625 +
                    im[((i  )*c+(j-1))]*.25 +
                    im[((i  )*c+(j  ))]*.375 +
                    im[((i  )*c+(j+1))]*.3125;
            COST_INC_ADD(3);
            COST_INC_MUL(4);
        }else{
            j = c-1; // [ ... -2 -1 0] 0 -1
            tmp_halfsize[(i*c2+(j/2))] =
                    im[((i  )*c+(j-2))]*.0625 +
                    im[((i  )*c+(j-1))]*.3125 +
                    im[((i  )*c+(j  ))]*.625;
            COST_INC_ADD(2);
            COST_INC_MUL(3);
        }
    }
    //vertical filter
    c = c2;
    for(int j = 0; j < c; j++){ //every column in tmp_halfsize = every 2nd column in im
        for(int i = 2; i < r-2; i+=2){ //every 2nd row in tmp_halfsize
            int i2 = i/2+((i-2) % 2); //TODO remove
            dst[(i2*c+j)] =
                    tmp_halfsize[((i-2)*c+(j  ))]*.0625 +
                    tmp_halfsize[((i-1)*c+(j  ))]*.25 +
                    tmp_halfsize[((i  )*c+(j  ))]*.375 +
                    tmp_halfsize[((i+1)*c+(j  ))]*.25 +
                    tmp_halfsize[((i+2)*c+(j  ))]*.0625;
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        //top edge
        int i = 0; // 1 0 [0 1 2 ... ]
        int i2 = i/2+((i-2) % 2);  //TODO remove
        dst[(i2*c+j)] =
                tmp_halfsize[((i+1)*c+(j  ))]*.3125 +
                tmp_halfsize[((i  )*c+(j  ))]*.625 +
                tmp_halfsize[((i+2)*c+(j  ))]*.0625;
        COST_INC_ADD(2);
        COST_INC_MUL(3);
        //bottom edge
        if((r-2) % 2 == 0){
            i = r-2; // [ ... -2 -1 0 1] 1
            i2 = i/2+((i-2) % 2);  //TODO remove
            dst[(i2*c+j)] =
                    tmp_halfsize[((i-2)*c+(j  ))]*.0625 +
                    tmp_halfsize[((i-1)*c+(j  ))]*.25 +
                    tmp_halfsize[((i  )*c+(j  ))]*.375 +
                    tmp_halfsize[((i+1)*c+(j  ))]*.3125;
            COST_INC_ADD(3);
            COST_INC_MUL(4);
        }else{
            i = r-1; // [ ... -2 -1 0] 0 -1
            i2 = i/2+((i-2) % 2);  //TODO remove
            dst[(i2*c+j)] =
                    tmp_halfsize[((i-2)*c+(j  ))]*.0625 +
                    tmp_halfsize[((i-1)*c+(j  ))]*.3125 +
                    tmp_halfsize[((i  )*c+(j  ))]*.625;
            COST_INC_ADD(2);
            COST_INC_MUL(3);
        }
    }
}

void upsample(double *im, uint32_t r, uint32_t c, uint32_t up_r, uint32_t up_c, double *tmp_halfsize, double *dst){
    // x
    for(int i = 0; i < r; i++){ //every 2nd line
        int j_half = 1;
        for(int j = 2; j < up_c-2; j++){
            tmp_halfsize[(i*up_c+j)] =
                    0.25*im[((i  )*c+(j_half-1))] +
                    1.5*im[((i )*c+(j_half  ))] +
                    0.25*im[((i  )*c+(j_half+1))];
            COST_INC_ADD(2);
            COST_INC_MUL(3);
            j++;
            tmp_halfsize[(i*up_c+j)] =
                    im[((i  )*c+(j_half))] +
                    im[((i  )*c+(j_half+1))];
            COST_INC_ADD(1);
            j_half++;
        }
        //left edge
        int j = 0; // 0 0 [0 1 2 ... ]
        tmp_halfsize[(i*up_c+j)] =
                1.75*im[(i  )*c] +
                0.25*im[(i  )*c+1];
        COST_INC_ADD(1);
        COST_INC_MUL(2);
        j = 1; // -1 [-1 0 1 2 ... ]
        tmp_halfsize[(i*up_c+j)] =
                im[((i  )*c+(0))] +
                im[((i  )*c+(1))];
        COST_INC_ADD(1);
        if(up_c % 2 == 0){
            //right edge
            j = up_c-2; // [ ... -2 -1 0 1] 0
            tmp_halfsize[(i*up_c+j)] =
                    0.25*im[((i  )*c+(j/2-1))] +
                    1.75*im[((i  )*c+(j/2  ))];
            COST_INC_ADD(1);
            COST_INC_MUL(2);
            j = up_c-1; // [ ... -2 -1 0] 0 0
            tmp_halfsize[(i*up_c+j)] =
                    2.0*im[((i  )*c+(j/2))];
            COST_INC_MUL(1);
        } else {
            //right edge (remaining)
            j = up_c-1; // [ ... -2 -1 0] 0 0
            tmp_halfsize[(i*up_c+j)] =
                    0.25*im[((i  )*c+(j/2-1))] +
                    1.75*im[((i  )*c+(j/2  ))];
            COST_INC_ADD(1);
            COST_INC_MUL(2);
        }
    }

    // y
    for(int j = 0; j < up_c; j++){ //all columns
        for(int i = 1; i < r-1; i++){
            dst[(2*i*up_c+j)] =
                    tmp_halfsize[((i-1)*up_c+(j  ))]*.0625 +
                    tmp_halfsize[((i  )*up_c+(j  ))]*.375 +
                    tmp_halfsize[((i+1)*up_c+(j  ))]*.0625;
            COST_INC_ADD(2);
            COST_INC_MUL(3);
            dst[((2*i+1)*up_c+j)] =
                    tmp_halfsize[((i)*up_c+(j  ))]*.25 +
                    tmp_halfsize[((i+1)*up_c+(j  ))]*.25;
            COST_INC_ADD(1);
            COST_INC_MUL(2);
        }
        //top edge
        int i = 0; // 0 0 [0 1 2 ... ]
        dst[(2*0*up_c+j)] =
                tmp_halfsize[((i  )*up_c+(j  ))]*.4375 +
                tmp_halfsize[((i+1)*up_c+(j  ))]*.0625;
        COST_INC_ADD(1);
        COST_INC_MUL(2);
        dst[((2*0+1)*up_c+j)] =
                tmp_halfsize[((i)*up_c+(j  ))]*.25 +
                tmp_halfsize[((i+1)*up_c+(j  ))]*.25;
        COST_INC_ADD(1);
        COST_INC_MUL(2);
        if(up_r % 2 == 0){
            //bottom edge
            i = r-1; // [ ... -2 -1 0 1] 0
            dst[((2*i)*up_c+j)] =
                    tmp_halfsize[((i-1)*up_c+(j  ))]*.0625 +
                    tmp_halfsize[((i  )*up_c+(j  ))]*.4375;
            COST_INC_ADD(1);
            COST_INC_MUL(2);
            // [ ... -2 -1 0] 0 0
            dst[((2*i+1)*up_c+j)] =
                    tmp_halfsize[((i)*up_c+(j  ))]*.5;
            COST_INC_MUL(1);
        }else{
            //bottom edge
            i = r-1; // [ ... -2 -1 0] 0 0
            dst[((2*i)*up_c+j)] =
                    tmp_halfsize[((i-1)*up_c+(j  ))]*.0625 +
                    tmp_halfsize[((i  )*up_c+(j  ))]*.4375;
            COST_INC_ADD(1);
            COST_INC_MUL(2);
        }
    }
}
