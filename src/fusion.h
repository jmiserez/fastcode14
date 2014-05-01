#ifndef FUSION_H
#define FUSION_H

#define FUSION_ALLOC_SUCCESS 0
#define FUSION_ALLOC_FAILURE -1

/**
 * @brief alloc_fusion allocates memory for internal data structures used by the fusion algorithm.
 * @param r height of the images
 * @param c width of the images
 * @param N number of images
 * @param segments object that holds pointers to the internally used memory segments
 * @return returns a number of >=0 if successful. <0 on error. In the latter case, you might want to check errno.
 */
int alloc_fusion(void** _segments, int w, int h, int N);

/**
 * @brief exposure_fusion
 * @param I
 * @param r
 * @param c
 * @param N
 * @param m
 * @param R
 * @param segments
 */
double* exposure_fusion(double** I, int w, int h, int N,
                        double contrast_parm, double sat_parm, double wexp_parm,
                        void* _segments);

void free_fusion( void* _segments, int N );

#endif // FUSION_H
