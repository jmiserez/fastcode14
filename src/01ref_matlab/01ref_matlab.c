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

//#define PRINTPYRAMIDS

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//TODO: rename these variables, preferably automatically
//TODO: remove length variables (W_len, etc.) here, we should not need them
//TODO: remove nlev
typedef struct {
    //W[n] is a weight map (1 value/pixel), there is one for each of the N images

    int w,h,N;
    size_t R_len; //TODO remove
    size_t W_len; //TODO remove
    size_t W_len2; //TODO remove
    size_t mono_len; //TODO remove
    size_t C_len; //TODO remove
    uint32_t nlev; //TODO remove (could be tricky)
    size_t Z_len; //TODO remove
    size_t S_len; //TODO remove
    size_t T_len; //TODO remove
    size_t Q_len; //TODO remove
    size_t U_len; //TODO remove
    size_t V_len; //TODO remove
    size_t pyrDisp_len; //TODO remove

    double *R;
    double **W;

    double *mono;
    //C is used as a temporary variable (1 value/pixel)
    double *C;

    //size fo pyramids

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

    //TODO: scratch space, reduce this as much as possible by inlining all operations

    //regular size scratch spaces
    double *Z;
    double *S;
    double *T;

    //large (+4 border) scratch spaces
    double *Q;
    double *U;
    double *V;

    //memory to print pyramids (solely used for debugging/verification)
    double *pyrDisp;

    //TODO: use these new names (or similar)
//    double** weight_matrices;
//    double* mono_matrix;
//    double* temp_matrix;
//    double* result;
} segments_t;

void exposure_fusion(double** I, int r, int c, int N, double m[3], double* R);

void contrast(double *im, uint32_t r, uint32_t c, double *mono, double *C);
void saturation(double *im, uint32_t npixels, double *C);
void well_exposedness(double *im, uint32_t npixels, double *C);
void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *Z, size_t Z_len, double *S, size_t S_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *Z, size_t Z_len, double *S, size_t S_len, double *T, size_t T_len, double *Q, size_t Q_len, double *U, size_t U_len, double *V, size_t V_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void reconstruct_laplacian_pyramid(uint32_t channels, uint32_t nlev, double *S, size_t S_len, double *Q, size_t Q_len, double *U, size_t U_len, double *V, size_t V_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst);
void downsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *Z, size_t Z_len, double *S, size_t S_len, uint32_t down_r, uint32_t down_c, double *dst);
void upsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *Q, size_t Q_len,  double *U, size_t U_len, double *V, size_t V_len, uint32_t up_r, uint32_t up_c, double *dst);
void display_pyramid(uint32_t channels, uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst);
void compact_display_pyramid(uint32_t channels, uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst);
void normalize_image(double *src, uint32_t channels, uint32_t r, uint32_t c, double *dst);

uint32_t compute_nlev(uint32_t r, uint32_t c);
int malloc_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr, uint32_t **pyr_r, uint32_t **pyr_c);
void free_pyramid(uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);

void rgb2gray(double *im, size_t npixels, double* dst);
void ones(double *dst, size_t len);
void ones_foreach(double **dst, size_t len, uint32_t N);
void zeros(double *dst, size_t len);
void zeros_foreach(double **dst, size_t len, uint32_t N);
void scalar_add(double *src, size_t src_len, double val, double *dst);
void scalar_mult(double *src, size_t src_len, double val, double *dst);
void scalar_div(double *src, size_t src_len, double val, double *dst);
void scalar_pow(double *src, size_t src_len, double val, double *dst);
void elementwise_mult(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_div(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_add(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_sub(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_sqrt(double *src, size_t src_len, double *dst);
void elementwise_copy(double *src, size_t src_len, double *dst);
void scalar_abs(double *src, size_t src_len, double *dst);
void conv3x3_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* fxy, double* dst);
void conv5x5separable_symmetric(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double *scratch, double* dst);
void conv5x5separable_replicate(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double *scratch, double* dst);

//TODO: move to separate (debug) c file
void store_image(char* path, double *R, uint32_t height, uint32_t width, uint32_t channels);

//
// Interface
//

int fusion_alloc(void** _segments, int w, int h, int N){

    segments_t *mem = malloc(sizeof(segments_t));
    if(mem == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->w = w;
    mem->h = h;
    mem->N = N;

    mem->R_len = w*h*3;
    mem->R = calloc(mem->R_len,sizeof(double));
    if(mem->R == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->W_len = N;
    mem->W_len2 = w*h;
    mem->W = calloc(mem->W_len,sizeof(double*));
    if(mem->W == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    for (int n = 0; n < N; n++){
        mem->W[n] = calloc(mem->W_len2,sizeof(double));
        if(mem->W[n] == NULL){
            return FUSION_ALLOC_FAILURE;
        }
    }

    mem->mono_len = w*h; //1 value/pixel
    mem->mono = calloc(mem->mono_len,sizeof(double));
    assert(mem->mono != NULL);

    mem->C_len = w*h;
    mem->C = calloc(mem->C_len,sizeof(double));
    if(mem->C == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->nlev = (uint32_t)(floor((log2(MIN(w,h)))));

    //TODO: allocate all pyramid memory as one chunk, not nlev chunks
    //TODO: remove writes to pyr_c and pyr_r!
    //      (as we are not allowed to already write to memory in this fusion_alloc() step).

    mem->pyrW = calloc(N,sizeof(double**));
    if(mem->pyrW == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrW_r = calloc(N,sizeof(double*));
    if(mem->pyrW_r == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrW_c = calloc(N,sizeof(double*));
    if(mem->pyrW_c == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->pyrI = calloc(N,sizeof(double**));
    if(mem->pyrI == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrI_r = calloc(N,sizeof(double*));
    if(mem->pyrI_r == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrI_c = calloc(N,sizeof(double*));
    if(mem->pyrI_c == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    if(malloc_pyramid(h,w,3,mem->nlev,&(mem->pyr), &(mem->pyr_r), &(mem->pyr_c)) != FUSION_ALLOC_SUCCESS){
        return FUSION_ALLOC_FAILURE;
    }
    for (int n = 0; n < N; n++){
        if(malloc_pyramid(h,w,1,mem->nlev,&(mem->pyrW[n]), &(mem->pyrW_r[n]), &(mem->pyrW_c[n])) != FUSION_ALLOC_SUCCESS){
            return FUSION_ALLOC_FAILURE;
        }
    }
    for (int n = 0; n < N; n++){
        if(malloc_pyramid(h,w,3,mem->nlev,&(mem->pyrI[n]), &(mem->pyrI_r[n]), &(mem->pyrI_c[n])) != FUSION_ALLOC_SUCCESS){
            return FUSION_ALLOC_FAILURE;
        }
    }

    //regular size scratch spaces
    mem->Z_len = w*h*3;
    mem->S_len = w*h*3;
    mem->T_len = w*h*3;

    mem->Z = calloc(mem->Z_len,sizeof(double));
    if(mem->Z == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->S = calloc(mem->S_len,sizeof(double));
    if(mem->S == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->T = calloc(mem->T_len,sizeof(double));
    if(mem->T == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    //large (+4 border) scratch spaces
    //(largest downsampled image + 1px border), upsampled 2x
    uint32_t largest_upsampled_r = (((h/2) + (h%2)) + 2) * 2;
    uint32_t largest_upsampled_c = (((w/2) + (w%2)) + 2) * 2;

    mem->Q_len = largest_upsampled_r*largest_upsampled_c*3;
    mem->U_len = largest_upsampled_r*largest_upsampled_c*3;
    mem->V_len = largest_upsampled_r*largest_upsampled_c*3;

    mem->Q = calloc(mem->Q_len,sizeof(double));
    if(mem->Q == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->U = calloc(mem->U_len,sizeof(double));
    if(mem->U == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->V = calloc(mem->V_len,sizeof(double));
    if(mem->V == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    mem->pyrDisp_len = (w*2)*h*3;
    mem->pyrDisp = calloc(mem->pyrDisp_len,sizeof(double));
    if(mem->pyrDisp == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    *_segments = mem;

    return FUSION_ALLOC_SUCCESS;
}

double* fusion_compute(double** I,
                        double contrast_parm, double sat_parm, double wexp_parm,
                        void* _segments){
    segments_t *mem = _segments;
    double* R = mem->R;

    int w = mem->w;
    int h = mem->h;
    int N = mem->N;

    int r = h;
    int c = w;

    size_t npixels = r*c;

    //W[n] is a weight map (1 value/pixel)
    //There is one for each of the N images
    size_t W_len2 = mem->W_len2;
    double** W = mem->W;
    assert(W != NULL);
    assert(mem->W_len == N);

#ifndef NDEBUG
    for (int n = 0; n < N; n++){
        assert(W[n] != NULL);
    }
#endif
PERF_FUNC_ENTER
    //C is used as a temporary variable (1 value/pixel)
    size_t C_len = mem->C_len;
    double* C = mem->C;
    assert(C != NULL);
    assert(W_len2 == C_len);

    double *mono = mem->mono;
    assert(mono != NULL);

    //for each image, calculate the weight maps
    for (int n = 0; n < N; n++){
        ones(W[n], W_len2);

        if(contrast_parm > 0){
            contrast(I[n],r,c,mono,C);
            scalar_pow(C,C_len,contrast_parm,C);
            elementwise_mult(W[n],W_len2,C,W[n]);
        }

        if(sat_parm > 0){
            saturation(I[n],npixels,C);
            scalar_pow(C,C_len,sat_parm,C);
            elementwise_mult(W[n],W_len2,C,W[n]);
        }

        if(wexp_parm > 0){
            well_exposedness(I[n],npixels,C);
            scalar_pow(C,C_len,wexp_parm,C);
            elementwise_mult(W[n],W_len2,C,W[n]);
        }

        scalar_add(W[n],W_len2,1.0E-12,W[n]);
    }

    //normalize weights: the total sum of weights for each pixel should be 1 across all N images
    elementwise_copy(W[0],W_len2,C);
    for (int n = 1; n < N; n++){
        elementwise_add(C,C_len,W[n],C);
    }
    for (int n = 0; n < N; n++){
        elementwise_div(W[n],W_len2,C,W[n]);
    }
PERF_FUNC_EXIT
#ifndef NDEBUG
    for(int i = 0; i < W_len2; i++){
        double sum_weight = 0;
        for (int n = 0; n < mem->W_len; n++){
            sum_weight += W[n][i];
        }
        assert(sum_weight >= 0.99 && sum_weight <= 1.01); //ensure all weights sum to one for each pixel
    }
#endif
#ifdef NO_PYRAMIDS
    return R;
#endif
#ifdef PRINTPYRAMIDS
    for(int i = 0; i < W_len; i++){
        char fname[50];
        sprintf(fname, "_W[%d].tif", i);
        store_image(fname, W[i], r, c, 1);
        printf("%s written to disk\n", fname);
    }
#endif

    uint32_t nlev = compute_nlev(r,c);
    assert(nlev != 0);

    //create empty pyramid
    double **pyr = mem->pyr;
    uint32_t *pyr_r = mem->pyr_r;
    uint32_t *pyr_c = mem->pyr_c;
    assert(pyr != NULL);
    assert(pyr_r != NULL);
    assert(pyr_c != NULL);

    for(int v = 0; v < nlev; v++){
        zeros(pyr[v],pyr_r[v]*pyr_c[v]*3);
    }

    //multiresolution blending
    double ***pyrW = mem->pyrW;
    uint32_t **pyrW_r = mem->pyrW_r;
    uint32_t **pyrW_c = mem->pyrW_c;
    assert(pyrW != NULL);
    assert(pyrW_r != NULL);
    assert(pyrW_c != NULL);

    double ***pyrI = mem->pyrI;
    uint32_t **pyrI_r = mem->pyrI_r;
    uint32_t **pyrI_c = mem->pyrI_c;
    assert(pyrI != NULL);
    assert(pyrI_r != NULL);
    assert(pyrI_c != NULL);

    //scratch space for gaussian/laplacian pyramid
    //TODO: optimize these away if possible
    size_t Z_len = mem->Z_len;
    double* Z = mem->Z;
    assert(Z != NULL);
    size_t S_len = mem->S_len;
    double* S = mem->S;
    assert(S != NULL);
    size_t T_len = mem->T_len;
    double* T = mem->T;
    assert(T != NULL);
    //how much scratch space is needed for upsampling?

    //(largest downsampled image + 1px border), upsampled 2x
//    uint32_t largest_upsampled_r = (((r/2) + (r%2)) + 2) * 2;
//    uint32_t largest_upsampled_c = (((c/2) + (c%2)) + 2) * 2;

    size_t Q_len = mem->Q_len;
    double* Q = mem->Q;
    assert(Q != NULL);
    size_t U_len = mem->U_len;
    double* U = mem->U;
    assert(U != NULL);
    size_t V_len = mem->V_len;
    double* V = mem->V;
    assert(V != NULL);

    //memory for display_pyramid (twice the width of the input images)
    size_t pyrDisp_len = mem->pyrDisp_len;
    double* pyrDisp = mem->pyrDisp;
    assert(pyrDisp != NULL);
    zeros(pyrDisp,pyrDisp_len);

    for (int n = 0; n < N; n++){
        //construct 1-channel gaussian pyramid from weights
        assert(pyrW[n] != NULL);
        assert(pyrW_r[n] != NULL);
        assert(pyrW_c[n] != NULL);
        gaussian_pyramid(W[n],r,c,1,nlev,Z,Z_len,S,S_len,pyrW[n],pyrW_r[n],pyrW_c[n]);

        //construct 3-channel laplacian pyramid from images
        assert(pyrI[n] != NULL);
        assert(pyrI_r[n] != NULL);
        assert(pyrI_c[n] != NULL);
        laplacian_pyramid(I[n],r,c,3,nlev,Z,Z_len,S,S_len,T,T_len,Q,Q_len,U,U_len,V,V_len,pyrI[n],pyrI_r[n],pyrI_c[n]);

        //weighted blend
        for(int v = 0; v < nlev; v++){
            for(int i = 0; i < pyrI_r[n][v]; i++){
                for(int j = 0; j < pyrI_c[n][v]; j++){
                    for(int k = 0; k < 3; k++){
                        pyr[v][(i*pyr_c[v]+j)*3+k] += pyrW[n][v][i*pyrI_c[n][v]+j] * pyrI[n][v][(i*pyrI_c[n][v]+j)*3+k];
                    }
                }
            }
        }
    }

    //reconstruct laplacian pyramid
    reconstruct_laplacian_pyramid(3,nlev,S,S_len,Q,Q_len,U,U_len,V,V_len,pyr,pyr_r,pyr_c,r,c,R);

#ifdef PRINTPYRAMIDS
    for(int i = 0; i < I_len; i++){
        char fname[50];

        //pyrW
        display_pyramid(1,nlev,pyrW[i],pyrW_r[i],pyrW_c[i],r,c*2,pyrDisp);
        sprintf(fname, "_pyrW[%d].tif", i);
        store_image(fname, pyrDisp, r, c*2, 1);
        printf("%s written to disk\n", fname);

        normalize_image(pyrDisp,1,r,c*2,pyrDisp);
        sprintf(fname, "_pyrW_norm[%d].tif", i);
        store_image(fname, pyrDisp, r, c*2, 1);
        printf("%s written to disk\n", fname);

        compact_display_pyramid(1,nlev,pyrW[i],pyrW_r[i],pyrW_c[i],r,c,pyrDisp);
        sprintf(fname, "_pyrW_compact[%d].tif", i);
        store_image(fname, pyrDisp, r, c, 1);
        printf("%s written to disk\n", fname);

        normalize_image(pyrDisp,1,r,c,pyrDisp);
        sprintf(fname, "_pyrW_compact_norm[%d].tif", i);
        store_image(fname, pyrDisp, r, c, 1);
        printf("%s written to disk\n", fname);

        //pyrI
        display_pyramid(3,nlev,pyrI[i],pyrI_r[i],pyrI_c[i],r,c*2,pyrDisp);
        sprintf(fname, "_pyrI[%d].tif", i);
        store_image(fname, pyrDisp, r, c*2, 3);
        printf("%s written to disk\n", fname);

        normalize_image(pyrDisp,3,r,c*2,pyrDisp);
        sprintf(fname, "_pyrI_norm[%d].tif", i);
        store_image(fname, pyrDisp, r, c*2, 3);
        printf("%s written to disk\n", fname);

        compact_display_pyramid(3,nlev,pyrI[i],pyrI_r[i],pyrI_c[i],r,c,pyrDisp);
        sprintf(fname, "_pyrI_compact[%d].tif", i);
        store_image(fname, pyrDisp, r, c, 3);
        printf("%s written to disk\n", fname);

        normalize_image(pyrDisp,3,r,c,pyrDisp);
        sprintf(fname, "_pyrI_compact_norm[%d].tif", i);
        store_image(fname, pyrDisp, r, c, 3);
        printf("%s written to disk\n", fname);
    }

    display_pyramid(3,nlev,pyr,pyr_r,pyr_c,r,c*2,pyrDisp);
    store_image("_pyr.tif", pyrDisp, r, c*2, 3);
    printf("_pyr.tif written to disk\n");
    normalize_image(pyrDisp,3,r,c*2,pyrDisp);
    store_image("_pyr_norm.tif", pyrDisp, r, c*2, 3);
    printf("_pyr_norm.tif written to disk\n");

    compact_display_pyramid(3,nlev,pyr,pyr_r,pyr_c,r,c,pyrDisp);
    store_image("_pyr_compact.tif", pyrDisp, r, c, 3);
    printf("_pyr_compact.tif written to disk\n");
    normalize_image(pyrDisp,3,r,c,pyrDisp);
    store_image("_pyr_compact_norm.tif", pyrDisp, r, c, 3);
    printf("_pyr_compact_norm.tif written to disk\n");
#endif

#ifndef NDEBUG
    printf("done\n");
#endif

    return R;
}

void fusion_free( void* _segments ){
    segments_t *mem = _segments;
    free(mem->mono);
    free(mem->C);
    free(mem->Z);
    free(mem->S);
    free(mem->T);
    free(mem->Q);
    free(mem->U);
    free(mem->V);

    free_pyramid(mem->nlev,mem->pyr,mem->pyr_r,mem->pyr_c);
    free(mem->pyrDisp);

    for (int n = 0; n < mem->W_len; n++){
        free(mem->W[n]);
        free_pyramid(mem->nlev,mem->pyrI[n],mem->pyrI_r[n],mem->pyrI_c[n]);
        free_pyramid(mem->nlev,mem->pyrW[n],mem->pyrW_r[n],mem->pyrW_c[n]);
    }
    free(mem->pyrI);
    free(mem->pyrI_r);
    free(mem->pyrI_c);
    free(mem->pyrW);
    free(mem->pyrW_r);
    free(mem->pyrW_c);
    free(mem->W);
    free(mem->R);
    free(_segments);
}

//
// Exposure Fusion functionality
//

void contrast(double *im, uint32_t r, uint32_t c, double *mono, double *C){
    //laplacian filter
    double h[] = {
        0.0, 1.0, 0.0,
        1.0, -4.0, 1.0,
        0.0, 1.0, 0.0
    };
    zeros(C, r*c);
    //for each image, calculate contrast measure on grayscale version of the image

    rgb2gray(im, r*c, mono);
    conv3x3_monochrome_replicate(mono,r,c,h,C);
}


void saturation(double *im, uint32_t npixels, double *C){
    //saturation is computed as the standard deviation of the color channels
    size_t C_len = npixels;
    zeros(C, C_len);

    // simple version
    for(int i = 0; i < npixels; i++){
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        double mu = (r + g + b) / 3.0;
        COST_INC_ADD(2);
        COST_INC_DIV(1);
        C[i] = sqrt((pow(r-mu,2) + pow(g-mu,2) + pow(b-mu,2))/3.0);
        COST_INC_SQRT(1);
        COST_INC_POW(3);
        COST_INC_DIV(1);
        COST_INC_ADD(5);
    }
}

void well_exposedness(double *im, uint32_t npixels, double *C){
    size_t C_len = npixels;
    zeros(C, C_len);

    double sig = 0.2;
    for(int i = 0; i < npixels; i++){
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        r = exp(-0.5*pow(r - 0.5,2) / pow(sig,2));
        g = exp(-0.5*pow(g - 0.5,2) / pow(sig,2));
        b = exp(-0.5*pow(b - 0.5,2) / pow(sig,2));
        COST_INC_ADD(3);
        COST_INC_MUL(3);
        COST_INC_DIV(3);
        COST_INC_POW(6);
        COST_INC_EXP(3);
        C[i] = r*g*b;
        COST_INC_MUL(2);
    }
}

void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *Z, size_t Z_len, double *S, size_t S_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    //pyr is an array of nlev arrays containing images of different (!) sizes
    //at this point pyr is already malloc-ed
    //pyr_r and pyr_c contain the sizes for each level

    assert(r == pyr_r[0]);
    assert(c == pyr_c[0]);

    //copy image to the finest level (note: MATLAB version is 1-indexed, here we use 0-indexing)

    elementwise_copy(im,r*c*channels,pyr[0]);

    size_t pyramid_filter_len = 5;
    double pyramid_filter[] = {.0625, .25, .375, .25, .0625};

    for(int v = 1; v < nlev; v++){
        //downsample image and store into level
        downsample(pyr[v-1],pyr_r[v-1],pyr_c[v-1],channels,pyramid_filter,pyramid_filter_len,Z,Z_len,S,S_len,pyr_r[v],pyr_c[v],pyr[v]);
    }
}

void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *Z, size_t Z_len, double *S, size_t S_len, double *T, size_t T_len, double *Q, size_t Q_len, double *U, size_t U_len, double *V, size_t V_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    //pyr is an array of nlev arrays containing images of different (!) sizes
    //at this point pyr is already malloc-ed
    //pyr_r and pyr_c contain the sizes for each level

    assert(r == pyr_r[0]);
    assert(c == pyr_c[0]);

    //copy image to the finest level (note: MATLAB version is 1-indexed, here we use 0-indexing)

    elementwise_copy(im,r*c*channels,pyr[0]);

    size_t pyramid_filter_len = 5;
    double pyramid_filter[] = {.0625, .25, .375, .25, .0625};

    //J = image
    elementwise_copy(im,r*c*channels,T); //TODO: optimize this copy away, can use pointer swaps

    uint32_t S_r = r;
    uint32_t S_c = c;
    uint32_t T_r = r;
    uint32_t T_c = c;
    assert(S_r*S_c*channels <= S_len);
    assert(T_r*T_c*channels <= T_len);

    for(int v = 0; v < nlev-1; v++){
        //downsample image T further, store in S
        S_r = pyr_r[v+1];
        S_c = pyr_c[v+1];
        downsample(T,T_r,T_c,channels,pyramid_filter,pyramid_filter_len,Z,Z_len,S,S_len,S_r,S_c,S);

        assert(T_r*T_c == pyr_r[v]*pyr_c[v]);
        //upsample image S, store temporarily in pyramid
        upsample(S,S_r,S_c,channels,pyramid_filter,pyramid_filter_len,Q,Q_len,U,U_len,V,V_len,pyr_r[v],pyr_c[v],pyr[v]);

        //subtract pyramid from T, store difference (T - upsampled image) in pyramid
        elementwise_sub(T,T_r*T_c*channels,pyr[v],pyr[v]);

        T_r = S_r;
        T_c = S_c;
        //continue with downsampled image remainder
        elementwise_copy(S,S_r*S_c*channels,T); //TODO: optimize this copy away, can use pointer swaps (handle first and last cases!)
    }
    //coarsest level, residual low pass imaggi
    assert(T_r*T_c*channels <= T_len);
    assert(T_r*T_c*channels == pyr_r[nlev-1]*pyr_c[nlev-1]*channels);
    elementwise_copy(T,T_r*T_c*channels,pyr[nlev-1]);
}

void reconstruct_laplacian_pyramid(uint32_t channels, uint32_t nlev, double *S, size_t S_len, double *Q, size_t Q_len, double *U, size_t U_len, double *V, size_t V_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst){

    size_t pyramid_filter_len = 5;
    double pyramid_filter[] = {.0625, .25, .375, .25, .0625};

    //copy low pass residual to dst
    elementwise_copy(pyr[nlev-1],pyr_r[nlev-1]*pyr_c[nlev-1]*channels,dst);

    for (int v = nlev-2; v >= 0; v--){
        //upsample to S
        upsample(dst,pyr_r[v+1],pyr_c[v+1],channels,pyramid_filter,pyramid_filter_len,Q,Q_len,U,U_len,V,V_len,pyr_r[v],pyr_c[v],S);

        //add current level to S, store in dst
        elementwise_add(pyr[v],pyr_r[v]*pyr_c[v]*channels,S,dst);
    }
}

void downsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *Z, size_t Z_len, double *S, size_t S_len, uint32_t down_r, uint32_t down_c, double *dst){
    assert(filter_len == 5);
    assert(filter != NULL);
    assert(S != NULL);
    assert(r*c*channels <= S_len);
    // [1] -> [1]
    // [1 2] -> [1]
    // [1 2 3] -> [1 3]
    // [1 2 3 4] -> [1 3]
    // width: W/2 + W%2
    assert(down_r == (r/2) + (r%2));
    assert(down_c == (c/2) + (c%2));

    //low pass filter
    //TODO: can optimize this, only need to calculate 1/4 of all pixels (can skip every second pixel as the result is never used)
    conv5x5separable_symmetric(im,r,c,channels,filter,filter,Z,S);
    //decimate, using every second entry
    for(int i = 0; i < down_r; i++){
        for(int j = 0; j < down_c; j++){
            for(int k = 0; k < channels; k++){
                assert(((i*2)*c+(j*2))*channels+k < r*c*channels); //bounds checking
                assert((i*down_c+j)*channels+k < down_r*down_c*channels); //bounds checking
                dst[(i*down_c+j)*channels+k] = S[((i*2)*c+(j*2))*channels+k];
            }
        }
    }
}

void upsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *Q, size_t Q_len, double *U, size_t U_len, double *V, size_t V_len, uint32_t up_r, uint32_t up_c, double *dst){
    assert(filter_len == 5);
    assert(filter != NULL);

    uint32_t padding = 1;

    //sizes with added 1 px border and size increase of 2x
    uint32_t r_upsampled = (r+2*padding)*2;
    uint32_t c_upsampled = (c+2*padding)*2;
    assert(up_r == r_upsampled - 4*padding - (up_r % 2));
    assert(up_c == c_upsampled - 4*padding - (up_c % 2));

    assert(U_len >= r_upsampled*c_upsampled*channels);
    assert(V_len == U_len);

    zeros(U,U_len);

    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            for(int k = 0; k < channels; k++){
                U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[(i*c+j)*channels+k];
                COST_INC_MUL(1);
            }
        }
    }
    //top row
    int i = -1;
    for(int j = 0; j < c; j++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i+1)*c+(j  ))*channels+k];
            COST_INC_MUL(1);
        }
    }
    //bottom row
    i = r;
    for(int j = 0; j < c; j++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i-1)*c+(j  ))*channels+k];
            COST_INC_MUL(1);
        }
    }
    //left edge
    int j = -1;
    for(int i = 0; i < r; i++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i  )*c+(j+1))*channels+k];
            COST_INC_MUL(1);
        }
    }
    //right edge
    j = c;
    for(int i = 0; i < r; i++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i  )*c+(j-1))*channels+k];
            COST_INC_MUL(1);
        }
    }
    //corners
    for(int k = 0; k < channels; k++){
        i = -1;
        j = -1;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i+1)*c+(j+1))*channels+k];
        j = c;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i+1)*c+(j-1))*channels+k];
        i = r;
        j = -1;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i-1)*c+(j+1))*channels+k];
        j = c;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i-1)*c+(j-1))*channels+k];
        COST_INC_MUL(4);
    }

    //blur
    conv5x5separable_replicate(U, r_upsampled, c_upsampled, channels, filter, filter, Q, V);

    //remove the border and copy result
    for(int i = 0; i < up_r; i++){
        for(int j = 0; j < up_c; j++){
            for(int k = 0; k < channels; k++){
                dst[(i*up_c+j)*channels+k] = V[((i+2)*c_upsampled+(j+2))*channels+k];
            }
        }
    }
}

/**
 * @brief Follows the Matlab implementation. Not the same as used in the paper.
 */
void display_pyramid(uint32_t channels, uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst){
    assert(r == pyr_r[0]);
    assert(c == pyr_c[0] * 2);

#ifndef NDEBUG
    //check if the image is high enough to show the whole pyramid in full
    int pyr_c_cumulative_sum = 0;
    for(int v = 0; v < nlev; v++){
        pyr_c_cumulative_sum += pyr_c[v];
        assert(pyr_c_cumulative_sum <= c);
    }
#endif

    //set to zeros
    for(int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            for (int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] = 0;
            }
        }
    }

    //fill pixels
    uint32_t c_offset = 0;
    uint32_t r_offset = 0;
    for(int v = 0; v < nlev; v++){
        for(int i = 0; i < pyr_r[v]; i++){
            for (int j = 0; j < pyr_c[v]; j++){
                for (int k = 0; k < channels; k++){
                    dst[((i+r_offset)*c+(j+c_offset))*channels+k] = pyr[v][(i*pyr_c[v]+j)*channels+k];
                }
            }
        }
        c_offset += pyr_c[v]; //offset
    }
}

/**
 * @brief Mimics the images shown in paper.
 */
void compact_display_pyramid(uint32_t channels, uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst){
    assert(r == pyr_r[0]);
    assert(c == pyr_c[0]);

    //set to zeros
    for(int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            for (int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] = 0;
            }
        }
    }

    //fill pixels
    uint32_t c_offset = 0;
    uint32_t r_offset = 0;
    for(int v = 0; v < nlev; v++){
        r_offset = r-pyr_r[v];
        for(int i = 0; i < pyr_r[v]; i++){
            for (int j = 0; j < pyr_c[v]; j++){
                for (int k = 0; k < channels; k++){
                    dst[((i+r_offset)*c+(j+c_offset))*channels+k] = pyr[v][(i*pyr_c[v]+j)*channels+k];
                }
            }
        }
        if(v > 0){
            c_offset += pyr_c[v]; //add column offset, but not the first time
        }
    }
}

void normalize_image(double *src, uint32_t channels, uint32_t r, uint32_t c, double *dst){
    //handle negative values
    double min_value = DBL_MAX;
    double max_value = -DBL_MAX;
    for(int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            for (int k = 0; k < channels; k++){
                double value = src[(i*c+j)*channels+k];

                if(value < min_value){
                    min_value = value;
                }

                if(value > max_value){
                    max_value = value;
                }
            }
        }
    }
    assert(min_value < DBL_MAX);
    assert(max_value > -DBL_MAX);

    if(min_value < 1.0E-5){
        for(int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                for (int k = 0; k < channels; k++){
                    double value = src[(i*c+j)*channels+k];
                    dst[(i*c+j)*channels+k] = (value - min_value) / (max_value - min_value);
                    COST_INC_ADD(2);
                    COST_INC_DIV(1);
                }
            }
        }
    } else {
        for(int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                for (int k = 0; k < channels; k++){
                    dst[(i*c+j)*channels+k] = src[(i*c+j)*channels+k];
                }
            }
        }
    }
}

//
// Helper functions
//

/**
 * Compute the highest possible pyramid
 */
uint32_t compute_nlev(uint32_t r, uint32_t c){
    COST_INC_EXP(2);
    return (uint32_t)(floor((log2(MIN(r,c)))));
}

/**
 * Allocate memory for the gaussian/laplacian pyramid at *pyr
 */
int malloc_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr, uint32_t **pyr_r, uint32_t **pyr_c){
    size_t pyr_len = nlev;
    *pyr = (double**) calloc(pyr_len,sizeof(double*));
    assert(*pyr != NULL);
    *pyr_r = (uint32_t*) calloc(nlev,sizeof(uint32_t));
    assert(*pyr_r != NULL);
    *pyr_c = (uint32_t*) calloc(nlev,sizeof(uint32_t));
    assert(*pyr_c != NULL);

    uint32_t r_level = r;
    uint32_t c_level = c;

    for(int i = 0; i < nlev; i++){
        (*pyr_r)[i] = r_level; //store dimension r at level i
        (*pyr_c)[i] = c_level; //store dimension c at level i

        size_t L_len = r_level*c_level*channels;
        double* L = calloc(L_len,sizeof(double));
        if(L == NULL){
            return FUSION_ALLOC_FAILURE;
        }

        (*pyr)[i] = L; //add entry to array of pointers to image levels

        // for next level, width if odd: (W-1)/2+1, otherwise: (W-1)/2
        r_level = r_level / 2 + (r_level % 2);
        c_level = c_level / 2 + (c_level % 2);
    }
    return FUSION_ALLOC_SUCCESS;
}

void free_pyramid(uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    for(int i = 0; i < nlev; i++){
        free(pyr[i]);
    }
    free(pyr_r);
    free(pyr_c);
    free(pyr);
}

//
// MATLAB-equivalent functionality
//

/**
 * @brief Implementation of the MATLAB rgb2gray function
 *
 * See: http://www.mathworks.com/help/images/ref/rgb2gray.html
 *
 * @param rgb Input image
 * @param npixels Size of image in pixels
 * @param gray (out) Output image
 */
void rgb2gray(double *im, size_t npixels, double* dst){
    for(int i = 0; i < npixels; i++){
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        dst[i] = 0.2989 * r + 0.5870 * g + 0.1140 * b; //retain luminance, discard hue and saturation
        COST_INC_ADD(2);
        COST_INC_MUL(3);
    }
}

void ones(double *dst, size_t len){
    for(int i = 0; i < len; i++){
        dst[i] = (double)1.0;
    }
}
void ones_foreach(double **dst, size_t len, uint32_t N){
    for(int i = 0; i < N; i++){
        ones(dst[i], len);
    }
}

void zeros(double *dst, size_t len){
    for(int i = 0; i < len; i++){
        dst[i] = (double)0.0;
    }
}
void zeros_foreach(double **dst, size_t len, uint32_t N){
    for(int i = 0; i < N; i++){
        zeros(dst[i],len);
    }
}

void scalar_add(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i] + val; COST_INC_ADD(1);
    }
}

void scalar_mult(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i] * val; COST_INC_MUL(1);
    }
}

void scalar_div(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i] / val; COST_INC_DIV(1);
    }
}

void scalar_pow(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = pow(fabs(src[i]),val); COST_INC_POW(1); COST_INC_ABS(1);
    }
}

void elementwise_mult(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] * src2[i]; COST_INC_MUL(1);
    }
}

void elementwise_div(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] / src2[i]; COST_INC_DIV(1);//beware of division by zero
    }
}

void elementwise_add(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] + src2[i]; COST_INC_ADD(1);
    }
}

void elementwise_copy(double *src, size_t src_len, double *dst){
//    memcpy(dst,src,src_len*sizeof(double));
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i];
    }
}

void elementwise_sub(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] - src2[i]; COST_INC_ADD(1);
    }
}

void elementwise_sqrt(double *src, size_t src_len, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = sqrt(src[i]); COST_INC_SQRT(1);
    }
}

void scalar_abs(double *src, size_t src_len, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = fabs(src[i]); COST_INC_ABS(1);
    }
}

/**
 * @brief convolution of a monochrome image with a 3x3 filter and border mode "replication"
 */
void conv3x3_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* f, double* dst){
    for(int i = 1; i < r-1; i++){
        for(int j = 1; j < c-1; j++){
            dst[i*c+j] =
                    im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[2] +
                    im[(i)  *c+(j-1)]*f[3] + im[(i)  *c+(j)]*f[4] + im[(i)  *c+(j+1)]*f[5] +
                    im[(i+1)*c+(j-1)]*f[6] + im[(i+1)*c+(j)]*f[7] + im[(i+1)*c+(j+1)]*f[8];
             COST_INC_ADD(6);
             COST_INC_MUL(8);
        }
    }
    //edges
    for(int i = 1; i < r-1; i++){
        int j = 0;
        dst[i*c+j] =
                im[(i-1)*c+(j)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[2] +
                im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j)]*f[4] + im[(i)  *c+(j+1)]*f[5] +
                im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j)]*f[7] + im[(i+1)*c+(j+1)]*f[8];
        COST_INC_ADD(6);
        COST_INC_MUL(8);
        j = c-1;
        dst[i*c+j] =
                im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j)]*f[2] +
                im[(i)  *c+(j-1)]*f[3] + im[(i)  *c+(j)]*f[4] + im[(i)  *c+(j)]*f[5] +
                im[(i+1)*c+(j-1)]*f[6] + im[(i+1)*c+(j)]*f[7] + im[(i+1)*c+(j)]*f[8];
        COST_INC_ADD(6);
        COST_INC_MUL(8);
    }
    for(int j = 1; j < c-1; j++){
        int i = 0;
        dst[i*c+j] =
                im[(i)  *c+(j-1)]*f[0] + im[(i)  *c+(j)]*f[1] + im[(i)  *c+(j+1)]*f[2] +
                im[(i)  *c+(j-1)]*f[3] + im[(i)  *c+(j)]*f[4] + im[(i)  *c+(j+1)]*f[5] +
                im[(i+1)*c+(j-1)]*f[6] + im[(i+1)*c+(j)]*f[7] + im[(i+1)*c+(j+1)]*f[8];
        COST_INC_ADD(6);
        COST_INC_MUL(8);
        i = r-1;
        dst[i*c+j] =
                im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[2] +
                im[(i)  *c+(j-1)]*f[3] + im[(i)  *c+(j)]*f[4] + im[(i)  *c+(j+1)]*f[5] +
                im[(i)  *c+(j-1)]*f[6] + im[(i)  *c+(j)]*f[7] + im[(i)  *c+(j+1)]*f[8];
        COST_INC_ADD(6);
        COST_INC_MUL(8);
    }
    //corners
    //top left
    int i = 0;
    int j = 0;
    dst[i*c+j] =
            im[(i  )*c+(j)]*f[0] + im[(i  )*c+(j)]*f[1] + im[(i  )*c+(j+1)]*f[2] +
            im[(i  )*c+(j)]*f[3] + im[(i  )*c+(j)]*f[4] + im[(i  )*c+(j+1)]*f[5] +
            im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j)]*f[7] + im[(i+1)*c+(j+1)]*f[8];
    COST_INC_ADD(6);
    COST_INC_MUL(8);
    //top right
    i = 0;
    j = c-1;
    dst[i*c+j] =
            im[(i  )*c+(j-1)]*f[0] + im[(i  )*c+(j)]*f[1] + im[(i  )*c+(j  )]*f[2] +
            im[(i  )*c+(j-1)]*f[3] + im[(i  )*c+(j)]*f[4] + im[(i  )*c+(j  )]*f[5] +
            im[(i+1)*c+(j-1)]*f[6] + im[(i+1)*c+(j)]*f[7] + im[(i+1)*c+(j  )]*f[8];
    COST_INC_ADD(6);
    COST_INC_MUL(8);
    //bottom left
    i = r-1;
    j = 0;
    dst[i*c+j] =
            im[(i-1)*c+(j  )]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[2] +
            im[(i  )*c+(j  )]*f[3] + im[(i  )*c+(j)]*f[4] + im[(i  )*c+(j+1)]*f[5] +
            im[(i  )*c+(j  )]*f[6] + im[(i  )*c+(j)]*f[7] + im[(i  )*c+(j+1)]*f[8];
    COST_INC_ADD(6);
    COST_INC_MUL(8);
    //bottom right
    i = r-1;
    j = c-1;
    dst[i*c+j] =
            im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j  )]*f[2] +
            im[(i  )*c+(j-1)]*f[3] + im[(i  )*c+(j)]*f[4] + im[(i  )*c+(j  )]*f[5] +
            im[(i  )*c+(j-1)]*f[6] + im[(i  )*c+(j)]*f[7] + im[(i  )*c+(j  )]*f[8];
    COST_INC_ADD(6);
    COST_INC_MUL(8);

}

/**
 * @brief convolution of a multi-channel image with a separable 5x5 filter and border mode "symmetric"
 */
void conv5x5separable_symmetric(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double* scratch, double* dst){
    //r is height (vertical), c is width (horizontal)

    //horizontal filter
    for(int i = 0; i < r; i++){ //all lines
        for(int j = 2; j < c-2; j++){
            for(int k = 0; k < channels; k++){
                scratch[(i*c+j)*channels+k] =
                        im[((i  )*c+(j-2))*channels+k]*fx[0] +
                        im[((i  )*c+(j-1))*channels+k]*fx[1] +
                        im[((i  )*c+(j  ))*channels+k]*fx[2] +
                        im[((i  )*c+(j+1))*channels+k]*fx[3] +
                        im[((i  )*c+(j+2))*channels+k]*fx[4];
                COST_INC_ADD(4);
                COST_INC_MUL(5);
            }
        }
        //left edge
        int j = 0; // 1 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j+1))*channels+k]*fx[0] +
                    im[((i  )*c+(j  ))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        j = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-1))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        //right edge
        j = c-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+1))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        j = c-1; // [ ... -2 -1 0] 0 -1
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j  ))*channels+k]*fx[3] +
                    im[((i  )*c+(j-1))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
    }
    //vertical filter
    for(int j = 0; j < c; j++){ //all columns
        for(int i = 2; i < r-2; i++){
            for(int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] =
                        scratch[((i-2)*c+(j  ))*channels+k]*fy[0] +
                        scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                        scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                        scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                        scratch[((i+2)*c+(j  ))*channels+k]*fy[4];
                COST_INC_ADD(4);
                COST_INC_MUL(5);
            }
        }
        //top edge
        int i = 0; // 1 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i+2)*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        i = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i+2)*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        //bottom edge
        i = r-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        i = r-1; // [ ... -2 -1 0] 0 -1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
    }
}

/**
 * @brief convolution of a multi-channel image with a separable 5x5 filter and border mode "replicate"
 */
void conv5x5separable_replicate(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double* scratch, double* dst){
    //r is height (vertical), c is width (horizontal)

    //horizontal filter
    for(int i = 0; i < r; i++){ //all lines
        for(int j = 2; j < c-2; j++){
            for(int k = 0; k < channels; k++){
                scratch[(i*c+j)*channels+k] =
                        im[((i  )*c+(j-2))*channels+k]*fx[0] +
                        im[((i  )*c+(j-1))*channels+k]*fx[1] +
                        im[((i  )*c+(j  ))*channels+k]*fx[2] +
                        im[((i  )*c+(j+1))*channels+k]*fx[3] +
                        im[((i  )*c+(j+2))*channels+k]*fx[4];
                COST_INC_ADD(4);
                COST_INC_MUL(5);
            }
        }
        //left edge
        int j = 0; // 0 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j  ))*channels+k]*fx[0] +
                    im[((i  )*c+(j  ))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        j = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-1))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        //right edge
        j = c-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+1))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        j = c-1; // [ ... -2 -1 0] 0 0
        for(int k = 0; k < channels; k++){
            scratch[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j  ))*channels+k]*fx[3] +
                    im[((i  )*c+(j  ))*channels+k]*fx[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
    }
    //vertical filter
    for(int j = 0; j < c; j++){ //all columns
        for(int i = 2; i < r-2; i++){
            for(int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] =
                        scratch[((i-2)*c+(j  ))*channels+k]*fy[0] +
                        scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                        scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                        scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                        scratch[((i+2)*c+(j  ))*channels+k]*fy[4];
                COST_INC_ADD(4);
                COST_INC_MUL(5);
            }
        }
        //top edge
        int i = 0; // 0 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i  )*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i+2)*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        i = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i+2)*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        //bottom edge
        i = r-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i+1)*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
        i = r-1; // [ ... -2 -1 0] 0 0
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    scratch[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    scratch[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[2] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[3] +
                    scratch[((i  )*c+(j  ))*channels+k]*fy[4];
            COST_INC_ADD(4);
            COST_INC_MUL(5);
        }
    }
}

//TODO: move to separate (debug) c file.
void store_image(char* path, double *R, uint32_t height, uint32_t width, uint32_t channels){
    assert(channels == 1 || channels == 3);
    uint32_t npixels = width*height;

    TIFF *out = TIFFOpen(path, "w");
    uint32_t* raster = (uint32_t*) _TIFFmalloc(npixels * sizeof(uint32_t));

    assert(raster != NULL);
    switch(channels){
        case 1:
            for(int i = 0; i < npixels; i++){
                uint32_t packed = 0;
                packed |= (uint32_t)(fmin(fmax(0.0, round(R[i]*255.0)), 255.0));
                packed |= ((uint32_t)(fmin(fmax(0.0, round(R[i]*255.0)), 255.0))) << 8;
                packed |= ((uint32_t)(fmin(fmax(0.0, round(R[i]*255.0)), 255.0))) << 16;
                packed |= ((uint32_t)255) << 24;
                raster[i] = packed;
            }
            break;
        case 3:
            for(int i = 0; i < npixels; i++){
                uint32_t packed = 0;
                packed |= (uint32_t)(fmin(fmax(0.0, round(R[i*3]*255.0)), 255.0));
                packed |= ((uint32_t)(fmin(fmax(0.0, round(R[i*3+1]*255.0)), 255.0))) << 8;
                packed |= ((uint32_t)(fmin(fmax(0.0, round(R[i*3+2]*255.0)), 255.0))) << 16;
                packed |= ((uint32_t)255) << 24;
                raster[i] = packed;
            }
            break;
    }

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, height);
    TIFFWriteEncodedStrip(out, 0, raster, width*height*sizeof(uint32_t));
    _TIFFfree(raster);
    TIFFClose(out);
}
