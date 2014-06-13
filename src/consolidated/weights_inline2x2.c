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
    COST_INC_POW(1);
    COST_INC_ABS(1);
    double t3 = re*ge*be;
    COST_INC_MUL(2);
    double t4 = fabs(t3);
    COST_INC_POW(1);
    COST_INC_ABS(1);

    double t5 = t2 * t4;
    COST_INC_MUL(1);

    double t6 = -4.0*grey+grey1+grey2+grey3+grey4;
    COST_INC_MUL(1);
    COST_INC_ADD(4);

    double t7 = fabs(t6);
    COST_INC_POW(1);
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

    int row_multiplier = 2;
    int rows_to_skip = (i_end - i_start) % row_multiplier;
    if(rows_to_skip != 0){
        //skip ahead a few rows
        for(int i = i_start; i < i_start+rows_to_skip; i++){
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
    i_start = i_start+rows_to_skip;
    j_start = j_start+cols_to_skip;

    double rW = 0.2989;
    double gW = 0.5870;
    double bW = 0.1140;

    //
    // Calculate 4 central values per iteration. Can reuse 4 values.
    //
    //      t2  t3 ---->
    //
    //  c1  c2  c3  c4  ---->
    //
    //  d1  d2  d3  d4  ---->
    //
    //      b2  b3  ----->
    //

    int i = i_start;
    for(i = i_start; i < i_end; i+=2){
        int j = j_start;

        int c1c = i*c+j-1; // == c3c of previous iteration
        int c2c = i*c+j; // == c4c of previous iteration
        int d1c = (i+1)*c+j-1; // == d3c of previous iteration
        int d2c = (i+1)*c+j; // == d4c of previous iteration

        int c1x3 = 3*c1c; // == c3c of previous iteration
        int c2x3 = 3*c2c; // == c4c of previous iteration
        int d1x3 = 3*d1c; // == d3c of previous iteration
        int d2x3 = 3*d2c; // == d4c of previous iteration

        double c3_prev =
                im[c1x3] * rW +
                im[c1x3+1] * gW +
                im[c1x3+2] * bW;
        double c4_prev_r = im[c2x3];
        double c4_prev_g = im[c2x3+1];
        double c4_prev_b = im[c2x3+2];
        double c4_prev =
                c4_prev_r * rW +
                c4_prev_g * gW +
                c4_prev_b * bW;
        double d3_prev =
                im[d1x3] * rW +
                im[d1x3+1] * gW +
                im[d1x3+2] * bW;
        double d4_prev_r = im[d2x3];
        double d4_prev_g = im[d2x3+1];
        double d4_prev_b = im[d2x3+2];
        double d4_prev =
                d4_prev_r * rW +
                d4_prev_g * gW +
                d4_prev_b * bW;
        for(j = j_start; j < j_end; j+=2){

            int t2c = (i-1)*c+j;
            int t3c = (i-1)*c+j+1;
            c2c = i*c+j; // == c4c of previous iteration
            int c3c = i*c+j+1;
            int c4c = i*c+j+2;
            d2c = (i+1)*c+j; // == d4c of previous iteration
            int d3c = (i+1)*c+j+1;
            int d4c = (i+1)*c+j+2;
            int b2c = (i+1)*c+j;
            int b3c = (i+1)*c+j+1;

            int t2x3 = 3*((i-1)*c+j);
            int t3x3 = 3*((i-1)*c+j+1);
            int c3x3 = 3*(i*c+j+1);
            int c4x3 = 3*(i*c+j+2);
            int d3x3 = 3*((i+1)*c+j+1);
            int d4x3 = 3*((i+1)*c+j+2);
            int b2x3 = 3*((i+2)*c+j);
            int b3x3 = 3*((i+2)*c+j+1);

            double t2_r = im[t2x3];
            double t2_g = im[t2x3+1];
            double t2_b = im[t2x3+2];
            double t3_r = im[t3x3];
            double t3_g = im[t3x3+1];
            double t3_b = im[t3x3+2];

            double c2_r = c4_prev_r;
            double c2_g = c4_prev_g;
            double c2_b = c4_prev_b;
            double c3_r = im[c3x3];
            double c3_g = im[c3x3+1];
            double c3_b = im[c3x3+2];
            double c4_r = im[c4x3];
            double c4_g = im[c4x3+1];
            double c4_b = im[c4x3+2];
            c4_prev_r = c4_r;
            c4_prev_g = c4_g;
            c4_prev_b = c4_b;

            double d2_r = d4_prev_r;
            double d2_g = d4_prev_g;
            double d2_b = d4_prev_b;
            double d3_r = im[d3x3];
            double d3_g = im[d3x3+1];
            double d3_b = im[d3x3+2];
            double d4_r = im[d4x3];
            double d4_g = im[d4x3+1];
            double d4_b = im[d4x3+2];
            d4_prev_r = d4_r;
            d4_prev_g = d4_g;
            d4_prev_b = d4_b;

            double b2_r = im[b2x3];
            double b2_g = im[b2x3+1];
            double b2_b = im[b2x3+2];
            double b3_r = im[b3x3];
            double b3_g = im[b3x3+1];
            double b3_b = im[b3x3+2];

            double t2 = rW * t2_r + gW * t2_g + bW * t2_b;
            double t3 = rW * t3_r + gW * t3_g + bW * t3_b;
            double c1 = c3_prev;
            double c2 = c4_prev;
            double c3 = rW * c3_r + gW * c3_g + bW * c3_b;
            double c4 = rW * c4_r + gW * c4_g + bW * c4_b;
            c3_prev = c3;
            c4_prev = c4;
            double d1 = d3_prev;
            double d2 = d4_prev;
            double d3 = rW * d3_r + gW * d3_g + bW * d3_b;
            double d4 = rW * d4_r + gW * d4_g + bW * d4_b;
            d3_prev = d3;
            d4_prev = d4;
            double b2 = rW * b2_r + gW * b2_g + bW * b2_b;
            double b3 = rW * b3_r + gW * b3_g + bW * b3_b;



            double c2_mu = (c2_r + c2_g + c2_b) / 3.0;
            double c3_mu = (c3_r + c3_g + c3_b) / 3.0;
            double d2_mu = (d2_r + d2_g + d2_b) / 3.0;
            double d3_mu = (d3_r + d3_g + d3_b) / 3.0;


            double c2_rmu = c2_r - c2_mu;
            double c2_gmu = c2_g - c2_mu;
            double c2_bmu = c2_b - c2_mu;
            double c3_rmu = c3_r - c3_mu;
            double c3_gmu = c3_g - c3_mu;
            double c3_bmu = c3_b - c3_mu;
            double d2_rmu = d2_r - d2_mu;
            double d2_gmu = d2_g - d2_mu;
            double d2_bmu = d2_b - d2_mu;
            double d3_rmu = d3_r - d3_mu;
            double d3_gmu = d3_g - d3_mu;
            double d3_bmu = d3_b - d3_mu;

            double c2_rmu2 = c2_rmu * c2_rmu;
            double c2_gmu2 = c2_gmu * c2_gmu;
            double c2_bmu2 = c2_bmu * c2_bmu;
            double c3_rmu2 = c3_rmu * c3_rmu;
            double c3_gmu2 = c3_gmu * c3_gmu;
            double c3_bmu2 = c3_bmu * c3_bmu;
            double d2_rmu2 = d2_rmu * d2_rmu;
            double d2_gmu2 = d2_gmu * d2_gmu;
            double d2_bmu2 = d2_bmu * d2_bmu;
            double d3_rmu2 = d3_rmu * d3_rmu;
            double d3_gmu2 = d3_gmu * d3_gmu;
            double d3_bmu2 = d3_bmu * d3_bmu;

            double c2_rz = c2_r - 0.5;
            double c2_gz = c2_g - 0.5;
            double c2_bz = c2_b - 0.5;
            double c3_rz = c3_r - 0.5;
            double c3_gz = c3_g - 0.5;
            double c3_bz = c3_b - 0.5;
            double d2_rz = d2_r - 0.5;
            double d2_gz = d2_g - 0.5;
            double d2_bz = d2_b - 0.5;
            double d3_rz = d3_r - 0.5;
            double d3_gz = d3_g - 0.5;
            double d3_bz = d3_b - 0.5;

            double c2_rz2 = c2_rz * c2_rz;
            double c2_gz2 = c2_gz * c2_gz;
            double c2_bz2 = c2_bz * c2_bz;
            double c3_rz2 = c3_rz * c3_rz;
            double c3_gz2 = c3_gz * c3_gz;
            double c3_bz2 = c3_bz * c3_bz;
            double d2_rz2 = d2_rz * d2_rz;
            double d2_gz2 = d2_gz * d2_gz;
            double d2_bz2 = d2_bz * d2_bz;
            double d3_rz2 = d3_rz * d3_rz;
            double d3_gz2 = d3_gz * d3_gz;
            double d3_bz2 = d3_bz * d3_bz;


            double c2_rexp = exp(-12.5*c2_rz2);
            double c2_gexp = exp(-12.5*c2_gz2);
            double c2_bexp = exp(-12.5*c2_bz2);
            double c3_rexp = exp(-12.5*c3_rz2);
            double c3_gexp = exp(-12.5*c3_gz2);
            double c3_bexp = exp(-12.5*c3_bz2);
            double d2_rexp = exp(-12.5*d2_rz2);
            double d2_gexp = exp(-12.5*d2_gz2);
            double d2_bexp = exp(-12.5*d2_bz2);
            double d3_rexp = exp(-12.5*d3_rz2);
            double d3_gexp = exp(-12.5*d3_gz2);
            double d3_bexp = exp(-12.5*d3_bz2);

            double c2_t1 = sqrt((c2_rmu2 + c2_gmu2 + c2_bmu2)/3.0);
            double c3_t1 = sqrt((c3_rmu2 + c3_gmu2 + c3_bmu2)/3.0);
            double d2_t1 = sqrt((d2_rmu2 + d2_gmu2 + d2_bmu2)/3.0);
            double d3_t1 = sqrt((d3_rmu2 + d3_gmu2 + d3_bmu2)/3.0);

            double c2_t2 = fabs(c2_t1);
            double c3_t2 = fabs(c3_t1);
            double d2_t2 = fabs(d2_t1);
            double d3_t2 = fabs(d3_t1);

            double c2_t3 = c2_rexp*c2_gexp*c2_bexp;
            double c3_t3 = c3_rexp*c3_gexp*c3_bexp;
            double d2_t3 = d2_rexp*d2_gexp*d2_bexp;
            double d3_t3 = d3_rexp*d3_gexp*d3_bexp;

            double c2_t4 = fabs(c2_t3);
            double c3_t4 = fabs(c3_t3);
            double d2_t4 = fabs(d2_t3);
            double d3_t4 = fabs(d3_t3);

            double c2_t5 = c2_t2 * c2_t4;
            double c3_t5 = c3_t2 * c3_t4;
            double d2_t5 = d2_t2 * d2_t4;
            double d3_t5 = d3_t2 * d3_t4;

            double c2_t6 = -4.0*c2 + c1 + t2 + c3 + d2;
            double c3_t6 = -4.0*c3 + c2 + t3 + c4 + d3;
            double d2_t6 = -4.0*d2 + d1 + c2 + d3 + b2;
            double d3_t6 = -4.0*d3 + d2 + c3 + d4 + b3;

            double c2_t7 = fabs(c2_t6);
            double c3_t7 = fabs(c3_t6);
            double d2_t7 = fabs(d2_t6);
            double d3_t7 = fabs(d3_t6);


            double c2_t8 = c2_t5 * c2_t7;
            double c3_t8 = c3_t5 * c3_t7;
            double d2_t8 = d2_t5 * d2_t7;
            double d3_t8 = d3_t5 * d3_t7;

            double c2_t9 = c2_t8 + 1.0E-12;
            double c3_t9 = c3_t8 + 1.0E-12;
            double d2_t9 = d2_t8 + 1.0E-12;
            double d3_t9 = d3_t8 + 1.0E-12;

            dst[c2c] = c2_t9;
            dst[c3c] = c3_t9;
            dst[d2c] = d2_t9;
            dst[d3c] = d3_t9;

            COST_INC_MUL(80);
            COST_INC_ADD(76);
            COST_INC_DIV(8);
            COST_INC_EXP(12);
            COST_INC_SQRT(4);
            COST_INC_ABS(12);
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
