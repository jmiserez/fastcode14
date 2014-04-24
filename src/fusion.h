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
int alloc_fusion(int r, int c, int N, void* segments);

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
void exposure_fusion(double** I, int r, int c, int N, double m[3], double* R, void* segments);

void free_fusion( fusion_segments_t* segments );

#endif // FUSION_H
