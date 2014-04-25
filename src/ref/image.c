#include "image.h"

/**
 * @brief Implementation of the MATLAB rgb2gray function
 *
 * See: http://www.mathworks.com/help/images/ref/rgb2gray.html
 *
 * @param dst (out) Output image
 * @param im Input image
 * @param npixels Size of image in pixels
 */
void rgb2gray( double* dst, double *im, size_t npixels ) {
    size_t i;
    for( i = 0; i < npixels; i++ ) {
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        //retain luminance, discard hue and saturation
        dst[i] = 0.2989 * r + 0.5870 * g + 0.1140 * b;
    }
}


