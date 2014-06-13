#ifndef BASICS_C
#define BASICS_C

#include "basics.h"

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

    mem->pyrI = calloc(mem->nimages,sizeof(double**));
    if(mem->pyrI == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrI_r = calloc(mem->nimages,sizeof(double*));
    if(mem->pyrI_r == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->pyrI_c = calloc(mem->nimages,sizeof(double*));
    if(mem->pyrI_c == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    if(malloc_pyramid(mem->r,mem->c,3,mem->nlev,&(mem->pyr), &(mem->pyr_r), &(mem->pyr_c),false) != FUSION_ALLOC_SUCCESS){
        return FUSION_ALLOC_FAILURE;
    }
    for (int n = 0; n < N; n++){
        if(malloc_pyramid(mem->r,mem->c,1,mem->nlev,&(mem->pyrW[n]), &(mem->pyrW_r[n]), &(mem->pyrW_c[n]),true) != FUSION_ALLOC_SUCCESS){
            return FUSION_ALLOC_FAILURE;
        }
    }
    for (int n = 0; n < N; n++){
        if(malloc_pyramid(mem->r,mem->c,3,mem->nlev,&(mem->pyrI[n]), &(mem->pyrI_r[n]), &(mem->pyrI_c[n]),false) != FUSION_ALLOC_SUCCESS){
            return FUSION_ALLOC_FAILURE;
        }
    }

    mem->tmp_fullsize = calloc(mem->max_upsampled_r * mem->max_upsampled_c * 3,sizeof(double));
    if(mem->tmp_fullsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp2_fullsize = calloc(mem->max_upsampled_r * mem->max_upsampled_c * 3,sizeof(double));
    if(mem->tmp2_fullsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp_halfsize = calloc(
                (mem->max_upsampled_r * mem->max_upsampled_c * 3)+1 / 2, sizeof(double));
    if(mem->tmp_halfsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp_quartsize = calloc(
                ((mem->max_upsampled_r * mem->max_upsampled_c * 3)+1 / 2)+1 / 2, sizeof(double));
    if(mem->tmp_quartsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }
    mem->tmp2_quartsize = calloc(
                ((mem->max_upsampled_r * mem->max_upsampled_c * 3)+1 / 2)+1 / 2, sizeof(double));
    if(mem->tmp2_quartsize == NULL){
        return FUSION_ALLOC_FAILURE;
    }

    *_segments = mem;

    return FUSION_ALLOC_SUCCESS;
}

/**
 * Allocate memory for the gaussian/laplacian pyramid at *pyr
 */
int malloc_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr, uint32_t **pyr_r, uint32_t **pyr_c, bool level0_is_ref){
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

        if(i != 0 || !level0_is_ref){
            size_t L_len = r_level*c_level*channels;
            double* L = calloc(L_len,sizeof(double));
            if(L == NULL){
                return FUSION_ALLOC_FAILURE;
            }
            (*pyr)[i] = L; //add entry to array of pointers to image levels
        } else {
            (*pyr)[i] = NULL;
        }
        // for next level, width if odd: (W-1)/2+1, otherwise: (W-1)/2
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

    free_pyramid(mem->nlev,mem->pyr,mem->pyr_r,mem->pyr_c,false);

    for (int n = 0; n < mem->nimages; n++){
        free(mem->W[n]);
        free_pyramid(mem->nlev,mem->pyrI[n],mem->pyrI_r[n],mem->pyrI_c[n], false);
        free_pyramid(mem->nlev,mem->pyrW[n],mem->pyrW_r[n],mem->pyrW_c[n], true);
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

void free_pyramid(uint32_t nlev, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, bool level0_is_ref){
    for(int i = 0; i < nlev; i++){
        if(i != 0 || !level0_is_ref){
            free(pyr[i]);
        }
    }
    free(pyr_r);
    free(pyr_c);
    free(pyr);
}

#endif
