#ifndef FUSION_H
#define FUSION_H

/**
  * fusion_segments_t holds the pointers to the memory segments internally used by this implements of the fusion algorithm.
  */
typedef struct {
} fusion_segments_t;

#define FUSION_ALLOC_SUCCESS 0
#define FUSION_ALLOC_FAILURE -1

/**
 * @brief alloc_fusion allocates memory for internal data structures used by the fusion algorithm.
 * @param r height of the images
 * @param c width of the images
 * @param N number of images
 * @param segments struct that holds pointers to the internally used memory segments
 * @return returns a number of >=0 if successful. <0 on error. In the latter case, you might want to check errno.
 */
int alloc_fusion(int r, int c, int N, fusion_segments_t *segments);

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
void exposure_fusion(double** I, int r, int c, int N, double m[3], double* R, fusion_segments_t* segments);

#endif // FUSION_H
