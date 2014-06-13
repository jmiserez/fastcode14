#ifndef FUSION_H
#define FUSION_H

#define FUSION_ALLOC_SUCCESS 0
#define FUSION_ALLOC_FAILURE -1

/**
 * @brief fusion_alloc allocates memory for internal data structures used by the 
 *        fusion algorithm.
 * @param r height of the images
 * @param c width of the images
 * @param N number of images
 * @param segments object that holds pointers to the internally used memory
 *        segments
 * @return returns a number of >=0 if successful. <0 on error. In the latter
 *         case, you might want to check errno.
 */
int fusion_alloc(void** _segments, int w, int h, int N);

/**
 * @brief exposure_fusion
 * @param I
 * @param m
 * @param R
 * @param segments
 */
double* fusion_compute(double** I, double contrast_parm, double sat_parm, double wexp_parm,
                        void* _segments);

void fusion_free( void* _segments );

#endif // FUSION_H
