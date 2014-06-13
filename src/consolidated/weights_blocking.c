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

FORCE_INLINE void convolve_block(double* im, int ii, int jj, int N, uint32_t r, uint32_t c, double* dst){

    //determine non-special case boundaries
    int i_start = ii;
    int i_end = ii+N;
    int j_start = jj;
    int j_end = jj+N;
    if( ii == 0){
        i_start = 1;
    }
    if( ii + N == r){
        i_end = r-1;
    }
    if( jj == 0){
        j_start = 1;
    }
    if( jj + N == c){
        j_end = c-1;
    }

    if(ii == 0){
        if( jj == 0){
            //top left
            int i = 0;
            int j = 0;
            int center =  i*c+j;
            int bottom = (i+1)*c+j;
            int right   = i*c+j+1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*center,3*center,3*right,3*bottom);
        }
        if( jj + N == c){
            //top right
            int i = 0;
            int j = c-1;
            int center =  i*c+j;
            int bottom = (i+1)*c+j;
            int left   = i*c+j-1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*center,3*left,3*center,3*bottom);
        }
        //need to fill in top row
        int i = 0;
        for(int j = j_start; j < j_end; j++){
            int center =  i*c+j;
            int bottom = (i+1)*c+j;
            int left   = i*c+j-1;
            int right   = i*c+j+1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*center,3*left,3*right,3*bottom);
        }
    }
    if(ii + N == r){
        if( jj == 0){
            //bottom left
            int i = r-1;
            int j = 0;
            int center =  i*c+j;
            int top    = (i-1)*c+j;
            int right   = i*c+j+1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*top,3*center,3*right,3*center);
        }
        if( jj + N == c){
            //bottom right
            int i = r-1;
            int j = c-1;
            int center =  i*c+j;
            int top    = (i-1)*c+j;
            int left   = i*c+j-1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*top,3*left,3*center,3*center);
        }
        //need to fill in bottom row
        int i = r-1;
        for(int j = j_start; j < j_end; j++){
            int center =  i*c+j;
            int top    = (i-1)*c+j;
            int left   = i*c+j-1;
            int right   = i*c+j+1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*top,3*left,3*right,3*center);
        }
    }
    if( jj == 0){
        //need to fill in left edge
        int j = 0;
        for(int i = i_start; i < i_end; i++){
            int center =  i*c+j;
            int top    = (i-1)*c+j;
            int bottom = (i+1)*c+j;
            int right   = i*c+j+1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*top,3*center,3*right,3*bottom);
        }
    }
    if( jj + N == c){
        //need to fill in right edge
        int j = c-1;
        for(int i = i_start; i < i_end; i++){
            int center =  i*c+j;
            int top    = (i-1)*c+j;
            int bottom = (i+1)*c+j;
            int left   = i*c+j-1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*top,3*left,3*center,3*bottom);
        }
    }

    //now we can safely process everything inside the _start and _end margins
    //read pixels from _start-1 until _end+1 and process all inside

    for(int i = i_start; i < i_end; i++){
        for(int j = j_start; j < j_end; j++){
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

}

void convolve_calculate(double* im, uint32_t r, uint32_t c, double* dst){
    int blocksize = NB;
    if(blocksize > r){
        blocksize = r;
    }
    if(blocksize > c){
        blocksize = c;
    }
    for(int ii = 0; ii < r; ii+=blocksize){
        for(int jj = 0; jj < c; jj+=blocksize){
            //NOTE: image size MUST be an exact multiple of blocksize
            convolve_block(im,ii,jj,blocksize,r,c,dst);
        }
    }
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
