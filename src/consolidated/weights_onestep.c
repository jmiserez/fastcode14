#ifndef WEIGHTS_C
#define WEIGHTS_C

#include "weights.h"

/**
  * Calculate all values in one step per pixel. Requires grabbing the neighboring pixels.
  */
FORCE_INLINE double single_pixel(double *im, int center, int top, int left, int right, int bottom){
    double r = im[center];
    double g = im[center+1];
    double b = im[center+2];

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
    COST_INC_ADD(10);
    COST_INC_MUL(15);
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
    double t2 = fabs(t1);
    COST_INC_ABS(1);
    double t3 = re*ge*be;
    COST_INC_MUL(2);
    double t4 = fabs(t3);
    COST_INC_ABS(1);

    double t5 = t2 * t4;
    COST_INC_MUL(1);

    double t6 = -4.0*grey+grey1+grey2+grey3+grey4;
    COST_INC_MUL(1);
    COST_INC_ADD(4);

    double t7 = fabs(t6);
    COST_INC_ABS(1);

    double t8 = t5 * t7;
    COST_INC_MUL(1);

    double t9 = t8 + 1.0E-12;
    COST_INC_ADD(1);

    return t9;
}

/**
 * 3x3 convolution with mode "replicate"
 */
void convolve_calculate(double* im, uint32_t r, uint32_t c, double* dst){
    for(int i = 1; i < r-1; i++){
        for(int j = 1; j < c-1; j++){
            int center =  i*c+j;
            int top    = (i-1)*c+j;
            int bottom = (i+1)*c+j;
            int left   = i*c+j-1;
            int right   = i*c+j+1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*top,3*left,3*right,3*bottom);
        }

    }
    //top row
    int i = 0;
    int j = 0;
    for(j = 1; j < c-1; j++){
        int center =  i*c+j;
        int bottom = (i+1)*c+j;
        int left   = i*c+j-1;
        int right   = i*c+j+1;
        dst[center] = single_pixel(im,
                                   3*center,
                                   3*center,3*left,3*right,3*bottom);
    }
    //bottom row
    i = r-1;
    for(j = 1; j < c-1; j++){
        int center =  i*c+j;
        int top    = (i-1)*c+j;
        int left   = i*c+j-1;
        int right   = i*c+j+1;
        dst[center] = single_pixel(im,
                                   3*center,
                                   3*top,3*left,3*right,3*center);
    }
    //left edge
    j = 0;
    for(i = 1; i < r-1; i++){
        int center =  i*c+j;
        int top    = (i-1)*c+j;
        int bottom = (i+1)*c+j;
        int right   = i*c+j+1;
        dst[center] = single_pixel(im,
                                   3*center,
                                   3*top,3*center,3*right,3*bottom);
    }
    //right edge
    j = c-1;
    for(i = 1; i < r-1; i++){
        int center =  i*c+j;
        int top    = (i-1)*c+j;
        int bottom = (i+1)*c+j;
        int left   = i*c+j-1;
        dst[center] = single_pixel(im,
                                   3*center,
                                   3*top,3*left,3*center,3*bottom);
    }
    //corners
    //top left
    i = 0;
    j = 0;
    int center =  i*c+j;
    int bottom = (i+1)*c+j;
    int right   = i*c+j+1;
    dst[center] = single_pixel(im,
                               3*center,
                               3*center,3*center,3*right,3*bottom);
    //top right
    i = 0;
    j = c-1;
    center =  i*c+j;
    bottom = (i+1)*c+j;
    int left   = i*c+j-1;
    dst[center] = single_pixel(im,
                               3*center,
                               3*center,3*left,3*center,3*bottom);
    //bottom left
    i = r-1;
    j = 0;
    center =  i*c+j;
    int top    = (i-1)*c+j;
    right   = i*c+j+1;
    dst[center] = single_pixel(im,
                               3*center,
                               3*top,3*center,3*right,3*center);
    //bottom right
    i = r-1;
    j = c-1;
    center =  i*c+j;
    top    = (i-1)*c+j;
    left   = i*c+j-1;
    dst[center] = single_pixel(im,
                               3*center,
                               3*top,3*left,3*center,3*center);
}

void weights(uint32_t nimages, uint32_t r, uint32_t c,
             double **I, double **W){
    //for each image, calculate the weight maps
    for (int n = 0; n < nimages; n++){
        convolve_calculate(I[n],r,c,W[n]);
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

#endif
