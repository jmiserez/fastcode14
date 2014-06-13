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

//    for(int i = i_start; i < i_end; i++){
//        for(int j = j_start; j < j_end; j++){
//            int center =  i*c+j;
//            int top    = (i-1)*c+j;
//            int bottom = (i+1)*c+j;
//            int left   = i*c+j-1;
//            int right   = i*c+j+1;
//            dst[center] = single_pixel(im,
//                                       3*center,
//                                       3*top,3*left,3*right,3*bottom);
//        }
//    }

//    int row_multiplier = 2;
//    int rows_to_skip = (i_end - i_start) % row_multiplier;
//    if(rows_to_skip != 0){
//        //skip ahead a few rows
//        for(int i = i_start; i < i_start+rows_to_skip; i++){
//            for(j = j_start; j < j_end; j++){
//                int center =  i*c+j;
//                int top    = (i-1)*c+j;
//                int bottom = (i+1)*c+j;
//                int left   = i*c+j-1;
//                int right   = i*c+j+1;
//                dst[center] = single_pixel(im,
//                                           3*center,
//                                           3*top,3*left,3*right,3*bottom);
//            }
//        }
//    }

    int col_multiplier = 2;
    int cols_to_skip = (j_end - j_start) % col_multiplier;
    if(cols_to_skip != 0){
        //skip ahead a few cols
        for(int i = i_start; i < i_end; i++){
            for(int j = j_start; j < j_start+cols_to_skip; j++){
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
    j_start = j_start+cols_to_skip;

    double rw = 0.2989;
    double gw = 0.5870;
    double bw = 0.1140;

    for(int i = i_start; i < i_end; i++){
        for(int j = j_start; j < j_end; j+=2){
            int center_at =  (i*c+j);
            int center =  3*center_at;
            int top    = 3*((i-1)*c+j);
            int bottom = 3*((i+1)*c+j);
            int left   = 3*(i*c+j-1);
            int right_at   = (i*c+j+1);
            int right   = 3*right_at;
            int center2 =  right;
            int center2_at =  right_at;
            int top2    = 3*((i-1)*c+j+1);
            int bottom2 = 3*((i+1)*c+j+1);
            int right2   = 3*(i*c+j+2);

            double r0 = im[center];
            double g0 = im[center+1];
            double b0 = im[center+2];

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

            double r02 = r3;
            double g02 = g3;
            double b02 = b3;
            double r12 = im[top2];
            double g12 = im[top2+1];
            double b12 = im[top2+2];

            double r32 = im[right2];
            double g32 = im[right2+1];
            double b32 = im[right2+2];
            double r42 = im[bottom2];
            double g42 = im[bottom2+1];
            double b42 = im[bottom2+2];

            double grey0 = rw * r0 + gw * g0 + bw * b0;
            double grey1 = rw * r1 + gw * g1 + bw * b1;
            double grey2 = rw * r2 + gw * g2 + bw * b2;
            double grey3 = rw * r3 + gw * g3 + bw * b3;
            double grey4 = rw * r4 + gw * g4 + bw * b4;

            double grey02 = grey3;
            double grey12 = rw * r12 + gw * g12 + bw * b12;
            double grey22 = grey0;
            double grey32 = rw * r32 + gw * g32 + bw * b32;
            double grey42 = rw * r42 + gw * g42 + bw * b42;
            COST_INC_ADD(16);
            COST_INC_MUL(24);
            double mu = (r0 + g0 + b0) / 3.0;
            double mu2 = (r02 + g02 + b02) / 3.0;
            COST_INC_ADD(4);
            COST_INC_DIV(2);
            double rmu = r0-mu;
            double gmu = g0-mu;
            double bmu = b0-mu;
            double rmu2 = r02-mu2;
            double gmu2 = g02-mu2;
            double bmu2 = b02-mu2;
            COST_INC_ADD(6);
            double rz = r0-0.5;
            double gz = g0-0.5;
            double bz = b0-0.5;
            double rz2 = r02-0.5;
            double gz2 = g02-0.5;
            double bz2 = b02-0.5;
            COST_INC_ADD(6);
            double rzrz = rz*rz;
            double gzgz = gz*gz;
            double bzbz = bz*bz;
            double rzrz2 = rz2*rz2;
            double gzgz2 = gz2*gz2;
            double bzbz2 = bz2*bz2;
            COST_INC_MUL(6);
            double re = exp(-12.5*rzrz);
            double ge = exp(-12.5*gzgz);
            double be = exp(-12.5*bzbz);
            double re2 = exp(-12.5*rzrz2);
            double ge2 = exp(-12.5*gzgz2);
            double be2 = exp(-12.5*bzbz2);
            COST_INC_EXP(6);
            COST_INC_MUL(6);
            double t1 = sqrt((rmu*rmu + gmu*gmu + bmu*bmu)/3.0);
            double t12 = sqrt((rmu2*rmu2 + gmu2*gmu2 + bmu2*bmu2)/3.0);
            COST_INC_SQRT(2);
            COST_INC_ADD(6);
            COST_INC_MUL(6);
            COST_INC_DIV(2);
            double t2 = fabs(t1);
            double t22 = fabs(t12);
            COST_INC_ABS(2);
            double t3 = re*ge*be;
            double t32 = re2*ge2*be2;
            COST_INC_MUL(4);
            double t4 = fabs(t3);
            double t42 = fabs(t32);
            COST_INC_ABS(2);

            double t5 = t2 * t4;
            double t52 = t22 * t42;
            COST_INC_MUL(2);

            double t6 = -4.0*grey0+grey1+grey2+grey3+grey4;
            double t62 = -4.0*grey02+grey12+grey22+grey32+grey42;
            COST_INC_MUL(2);
            COST_INC_ADD(8);

            double t7 = fabs(t6);
            double t72 = fabs(t62);
            COST_INC_ABS(2);

            double t8 = t5 * t7;
            double t82 = t52 * t72;
            COST_INC_MUL(2);

            double t9 = t8 + 1.0E-12;
            double t92 = t82 + 1.0E-12;
            COST_INC_ADD(2);

            dst[center_at] = t9;
            dst[center2_at] = t92;
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
