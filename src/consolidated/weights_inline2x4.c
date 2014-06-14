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
    COST_INC_LOAD(15);

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
//        printf("skipped %d rows, i_start=%d, i_end=%d\n",rows_to_skip,i_start,i_end);
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

    int col_multiplier = 4;
    int cols_to_skip = (j_end - j_start) % col_multiplier;
    if(cols_to_skip != 0){
        //skip ahead a few cols
//        printf("skipped %d cols, j_start=%d, j_end=%d\n",cols_to_skip,j_start,j_end);
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
    // Calculate 8 central values per iteration. Can reuse 4 values.
    //
    //      t2  t3  ---->
    //
    //  c1  c2  c3  c4  ---->
    //
    //  d1  d2  d3  d4  ---->
    //
    //      b2  b3  ---->
    //
    //
    //
    //
    //  t1  t2  t3  t4  ---->
    //      rgb rgb
    //
    //  c1  c2  c3  c4  ---->
    //  rgb rgb rgb rgb
    //
    //  d1  d2  d3  d4  ---->
    //  rgb rgb rgb rgb
    //
    //  b1  b2  b3  b4  ---->
    //      rgb rgb
    //
    //  load in 3 ops:
    //  111 122 223 333
    //
    //
    //      t2  t3  t4  t5---->
    //
    //  c1 |c2  c3  c4  c5| c6  ---->
    //
    //  d1 |d2  d3  d4  d5| d6  ---->
    //
    //      b2  b3  b4  b5---->
    //
    //
    //  Reuse:
    //
    //      t2  t3  t4  t5      ---->
    //                      t2  t3  t4  t5
    //
    //  c1  c2  c3  c4  c5  c6  ---->
    //                  c1 |c2  c3  c4  c5| c6
    //
    //  d1  d2  d3  d4  d5  d6  ---->
    //                  d1 |d2  d3  d4  d5| d6
    //
    //      b2  b3  b4  b5      ---->
    //                      b2  b3  b4  b5
    //
    //
    //  Load c3456, d3456, t2345, b2345 every iteration. Shuffle the rest.
    //

    int i = i_start;
    for(i = i_start; i < i_end; i+=2){
        int j = j_start;

        //preload c3456, d3456 for first iteration
        int c1_x3 = 3*(i*c+j-1);
        int d1_x3 = 3*((i+1)*c+j-1);

        //load c3456
        double c3456_0 = 0;
        double c3456_1 = 0;
        double c3456_2 = 0;
        double c3456_3 = 0;
        double c3456_4 = 0;
        double c3456_5 = 0;
        double c3456_6 = im[c1_x3];
        double c3456_7 = im[c1_x3+1];
        double c3456_8 = im[c1_x3+2];
        double c3456_9 = im[c1_x3+3];
        double c3456_10 = im[c1_x3+4];
        double c3456_11 = im[c1_x3+5];

        //load d3456
        double d3456_0 = 0;
        double d3456_1 = 0;
        double d3456_2 = 0;
        double d3456_3 = 0;
        double d3456_4 = 0;
        double d3456_5 = 0;
        double d3456_6 = im[d1_x3];
        double d3456_7 = im[d1_x3+1];
        double d3456_8 = im[d1_x3+2];
        double d3456_9 = im[d1_x3+3];
        double d3456_10 = im[d1_x3+4];
        double d3456_11 = im[d1_x3+5];

        double c3456_grey_0 = 0;
        double c3456_grey_1 = 0;
        double c3456_grey_2 =
                rW * c3456_6 +
                gW * c3456_7 +
                bW * c3456_8;
        double c3456_grey_3 =
                rW * c3456_9 +
                gW * c3456_10 +
                bW * c3456_11;

        double d3456_grey_0 = 0;
        double d3456_grey_1 = 0;
        double d3456_grey_2 =
                rW * d3456_6 +
                gW * d3456_7 +
                bW * d3456_8;
        double d3456_grey_3 =
                rW * d3456_9 +
                gW * d3456_10 +
                bW * d3456_11;
        for(j = j_start; j < j_end; j+=4){

            //start indices
            int c2c = i*c+j;
            int d2c = (i+1)*c+j;
            int t2_x3 = 3*((i-1)*c+j);
            int c3_x3 = 3*(i*c+j+1);
            int d3_x3 = 3*((i+1)*c+j+1);
            int b2_x3 = 3*((i+2)*c+j);

            //load t2345
            double t2345_0 = im[t2_x3];
            double t2345_1 = im[t2_x3+1];
            double t2345_2 = im[t2_x3+2];
            double t2345_3 = im[t2_x3+3];
            double t2345_4 = im[t2_x3+4];
            double t2345_5 = im[t2_x3+5];
            double t2345_6 = im[t2_x3+6];
            double t2345_7 = im[t2_x3+7];
            double t2345_8 = im[t2_x3+8];
            double t2345_9 = im[t2_x3+9];
            double t2345_10 = im[t2_x3+10];
            double t2345_11 = im[t2_x3+11];

            //load c12
//            double c12_0 = c3456_6;
//            double c12_1 = c3456_7;
//            double c12_2 = c3456_8;
            double c12_3 = c3456_9;
            double c12_4 = c3456_10;
            double c12_5 = c3456_11;

            //load d12
//            double d12_0 = d3456_6;
//            double d12_1 = d3456_7;
//            double d12_2 = d3456_8;
            double d12_3 = d3456_9;
            double d12_4 = d3456_10;
            double d12_5 = d3456_11;

            //load c3456
            c3456_0 = im[c3_x3];
            c3456_1 = im[c3_x3+1];
            c3456_2 = im[c3_x3+2];
            c3456_3 = im[c3_x3+3];
            c3456_4 = im[c3_x3+4];
            c3456_5 = im[c3_x3+5];
            c3456_6 = im[c3_x3+6];
            c3456_7 = im[c3_x3+7];
            c3456_8 = im[c3_x3+8];
            c3456_9 = im[c3_x3+9];
            c3456_10 = im[c3_x3+10];
            c3456_11 = im[c3_x3+11];

            //load d3456
            d3456_0 = im[d3_x3];
            d3456_1 = im[d3_x3+1];
            d3456_2 = im[d3_x3+2];
            d3456_3 = im[d3_x3+3];
            d3456_4 = im[d3_x3+4];
            d3456_5 = im[d3_x3+5];
            d3456_6 = im[d3_x3+6];
            d3456_7 = im[d3_x3+7];
            d3456_8 = im[d3_x3+8];
            d3456_9 = im[d3_x3+9];
            d3456_10 = im[d3_x3+10];
            d3456_11 = im[d3_x3+11];

            //load b2345
            double b2345_0 = im[b2_x3];
            double b2345_1 = im[b2_x3+1];
            double b2345_2 = im[b2_x3+2];
            double b2345_3 = im[b2_x3+3];
            double b2345_4 = im[b2_x3+4];
            double b2345_5 = im[b2_x3+5];
            double b2345_6 = im[b2_x3+6];
            double b2345_7 = im[b2_x3+7];
            double b2345_8 = im[b2_x3+8];
            double b2345_9 = im[b2_x3+9];
            double b2345_10 = im[b2_x3+10];
            double b2345_11 = im[b2_x3+11];

            //copy old greyscale values
            double c12_grey_0 = c3456_grey_2;
            double c12_grey_1 = c3456_grey_3;
            double d12_grey_0 = d3456_grey_2;
            double d12_grey_1 = d3456_grey_3;

            //calculate greyscales
            double t2345_grey_0 =
                    rW * t2345_0 +
                    gW * t2345_1 +
                    bW * t2345_2;
            double t2345_grey_1 =
                    rW * t2345_3 +
                    gW * t2345_4 +
                    bW * t2345_5;
            double t2345_grey_2 =
                    rW * t2345_6 +
                    gW * t2345_7 +
                    bW * t2345_8;
            double t2345_grey_3 =
                    rW * t2345_9 +
                    gW * t2345_10 +
                    bW * t2345_11;

            double b2345_grey_0 =
                    rW * b2345_0 +
                    gW * b2345_1 +
                    bW * b2345_2;
            double b2345_grey_1 =
                    rW * b2345_3 +
                    gW * b2345_4 +
                    bW * b2345_5;
            double b2345_grey_2 =
                    rW * b2345_6 +
                    gW * b2345_7 +
                    bW * b2345_8;
            double b2345_grey_3 =
                    rW * b2345_9 +
                    gW * b2345_10 +
                    bW * b2345_11;

            c3456_grey_0 =
                    rW * c3456_0 +
                    gW * c3456_1 +
                    bW * c3456_2;
            c3456_grey_1 =
                    rW * c3456_3 +
                    gW * c3456_4 +
                    bW * c3456_5;
            c3456_grey_2 =
                    rW * c3456_6 +
                    gW * c3456_7 +
                    bW * c3456_8;
            c3456_grey_3 =
                    rW * c3456_9 +
                    gW * c3456_10 +
                    bW * c3456_11;

            d3456_grey_0 =
                    rW * d3456_0 +
                    gW * d3456_1 +
                    bW * d3456_2;
            d3456_grey_1 =
                    rW * d3456_3 +
                    gW * d3456_4 +
                    bW * d3456_5;
            d3456_grey_2 =
                    rW * d3456_6 +
                    gW * d3456_7 +
                    bW * d3456_8;
            d3456_grey_3 =
                    rW * d3456_9 +
                    gW * d3456_10 +
                    bW * d3456_11;


            //shuffle for outputs
            double c2345_0 = c12_3;
            double c2345_1 = c12_4;
            double c2345_2 = c12_5;
            double c2345_3 = c3456_0;
            double c2345_4 = c3456_1;
            double c2345_5 = c3456_2;
            double c2345_6 = c3456_3;
            double c2345_7 = c3456_4;
            double c2345_8 = c3456_5;
            double c2345_9 = c3456_6;
            double c2345_10 = c3456_7;
            double c2345_11 = c3456_8;

            double d2345_0 = d12_3;
            double d2345_1 = d12_4;
            double d2345_2 = d12_5;
            double d2345_3 = d3456_0;
            double d2345_4 = d3456_1;
            double d2345_5 = d3456_2;
            double d2345_6 = d3456_3;
            double d2345_7 = d3456_4;
            double d2345_8 = d3456_5;
            double d2345_9 = d3456_6;
            double d2345_10 = d3456_7;
            double d2345_11 = d3456_8;

            //shuffle greyscale for outputs
            double c2345_grey_0 = c12_grey_1;
            double c2345_grey_1 = c3456_grey_0;
            double c2345_grey_2 = c3456_grey_1;
            double c2345_grey_3 = c3456_grey_2;

            double d2345_grey_0 = d12_grey_1;
            double d2345_grey_1 = d3456_grey_0;
            double d2345_grey_2 = d3456_grey_1;
            double d2345_grey_3 = d3456_grey_2;

            //do convolution, center
            double x1_c2345_0 = c2345_grey_0 * -4.0;
            double x1_c2345_1 = c2345_grey_1 * -4.0;
            double x1_c2345_2 = c2345_grey_2 * -4.0;
            double x1_c2345_3 = c2345_grey_3 * -4.0;

            double x1_d2345_0 = d2345_grey_0 * -4.0;
            double x1_d2345_1 = d2345_grey_1 * -4.0;
            double x1_d2345_2 = d2345_grey_2 * -4.0;
            double x1_d2345_3 = d2345_grey_3 * -4.0;

            //do convolution, left
            double x2_c2345_0 = x1_c2345_0 + c12_grey_0;
            double x2_c2345_1 = x1_c2345_1 + c12_grey_1;
            double x2_c2345_2 = x1_c2345_2 + c2345_grey_1;
            double x2_c2345_3 = x1_c2345_3 + c2345_grey_2;

            double x2_d2345_0 = x1_d2345_0 + d12_grey_0;
            double x2_d2345_1 = x1_d2345_1 + d12_grey_1;
            double x2_d2345_2 = x1_d2345_2 + d2345_grey_1;
            double x2_d2345_3 = x1_d2345_3 + d2345_grey_2;

            //do convolution, right
            double x3_c2345_0 = x2_c2345_0 + c3456_grey_0;
            double x3_c2345_1 = x2_c2345_1 + c3456_grey_1;
            double x3_c2345_2 = x2_c2345_2 + c3456_grey_2;
            double x3_c2345_3 = x2_c2345_3 + c3456_grey_3;

            double x3_d2345_0 = x2_d2345_0 + d3456_grey_0;
            double x3_d2345_1 = x2_d2345_1 + d3456_grey_1;
            double x3_d2345_2 = x2_d2345_2 + d3456_grey_2;
            double x3_d2345_3 = x2_d2345_3 + d3456_grey_3;

            //do convolution, top
            double x4_c2345_0 = x3_c2345_0 + t2345_grey_0;
            double x4_c2345_1 = x3_c2345_1 + t2345_grey_1;
            double x4_c2345_2 = x3_c2345_2 + t2345_grey_2;
            double x4_c2345_3 = x3_c2345_3 + t2345_grey_3;

            double x4_d2345_0 = x3_d2345_0 + c2345_grey_0;
            double x4_d2345_1 = x3_d2345_1 + c2345_grey_1;
            double x4_d2345_2 = x3_d2345_2 + c2345_grey_2;
            double x4_d2345_3 = x3_d2345_3 + c2345_grey_3;

            //do convolution, bottom
            double x5_c2345_0 = x4_c2345_0 + d2345_grey_0;
            double x5_c2345_1 = x4_c2345_1 + d2345_grey_1;
            double x5_c2345_2 = x4_c2345_2 + d2345_grey_2;
            double x5_c2345_3 = x4_c2345_3 + d2345_grey_3;

            double x5_d2345_0 = x4_d2345_0 + b2345_grey_0;
            double x5_d2345_1 = x4_d2345_1 + b2345_grey_1;
            double x5_d2345_2 = x4_d2345_2 + b2345_grey_2;
            double x5_d2345_3 = x4_d2345_3 + b2345_grey_3;

            double x6_c2345_0 = fabs(x5_c2345_0);
            double x6_c2345_1 = fabs(x5_c2345_1);
            double x6_c2345_2 = fabs(x5_c2345_2);
            double x6_c2345_3 = fabs(x5_c2345_3);

            double x6_d2345_0 = fabs(x5_d2345_0);
            double x6_d2345_1 = fabs(x5_d2345_1);
            double x6_d2345_2 = fabs(x5_d2345_2);
            double x6_d2345_3 = fabs(x5_d2345_3);

            //convolution done, result in x6

            // shuffle
            double r_c2345_0 = c2345_0;
            double r_c2345_1 = c2345_3;
            double r_c2345_2 = c2345_6;
            double r_c2345_3 = c2345_9;

            double g_c2345_0 = c2345_1;
            double g_c2345_1 = c2345_4;
            double g_c2345_2 = c2345_7;
            double g_c2345_3 = c2345_10;

            double b_c2345_0 = c2345_2;
            double b_c2345_1 = c2345_5;
            double b_c2345_2 = c2345_8;
            double b_c2345_3 = c2345_11;

            double r_d2345_0 = d2345_0;
            double r_d2345_1 = d2345_3;
            double r_d2345_2 = d2345_6;
            double r_d2345_3 = d2345_9;

            double g_d2345_0 = d2345_1;
            double g_d2345_1 = d2345_4;
            double g_d2345_2 = d2345_7;
            double g_d2345_3 = d2345_10;

            double b_d2345_0 = d2345_2;
            double b_d2345_1 = d2345_5;
            double b_d2345_2 = d2345_8;
            double b_d2345_3 = d2345_11;

            double rgbsum_c2345_0 = r_c2345_0 + g_c2345_0 + b_c2345_0;
            double rgbsum_c2345_1 = r_c2345_1 + g_c2345_1 + b_c2345_1;
            double rgbsum_c2345_2 = r_c2345_2 + g_c2345_2 + b_c2345_2;
            double rgbsum_c2345_3 = r_c2345_3 + g_c2345_3 + b_c2345_3;

            double rgbsum_d2345_0 = r_d2345_0 + g_d2345_0 + b_d2345_0;
            double rgbsum_d2345_1 = r_d2345_1 + g_d2345_1 + b_d2345_1;
            double rgbsum_d2345_2 = r_d2345_2 + g_d2345_2 + b_d2345_2;
            double rgbsum_d2345_3 = r_d2345_3 + g_d2345_3 + b_d2345_3;

            double mu_c2345_0 = rgbsum_c2345_0 / 3.0;
            double mu_c2345_1 = rgbsum_c2345_1 / 3.0;
            double mu_c2345_2 = rgbsum_c2345_2 / 3.0;
            double mu_c2345_3 = rgbsum_c2345_3 / 3.0;

            double mu_d2345_0 = rgbsum_d2345_0 / 3.0;
            double mu_d2345_1 = rgbsum_d2345_1 / 3.0;
            double mu_d2345_2 = rgbsum_d2345_2 / 3.0;
            double mu_d2345_3 = rgbsum_d2345_3 / 3.0;


            double rmu_c2345_0 = r_c2345_0 - mu_c2345_0;
            double rmu_c2345_1 = r_c2345_1 - mu_c2345_1;
            double rmu_c2345_2 = r_c2345_2 - mu_c2345_2;
            double rmu_c2345_3 = r_c2345_3 - mu_c2345_3;

            double gmu_c2345_0 = g_c2345_0 - mu_c2345_0;
            double gmu_c2345_1 = g_c2345_1 - mu_c2345_1;
            double gmu_c2345_2 = g_c2345_2 - mu_c2345_2;
            double gmu_c2345_3 = g_c2345_3 - mu_c2345_3;

            double bmu_c2345_0 = b_c2345_0 - mu_c2345_0;
            double bmu_c2345_1 = b_c2345_1 - mu_c2345_1;
            double bmu_c2345_2 = b_c2345_2 - mu_c2345_2;
            double bmu_c2345_3 = b_c2345_3 - mu_c2345_3;

            double rmu_d2345_0 = r_d2345_0 - mu_d2345_0;
            double rmu_d2345_1 = r_d2345_1 - mu_d2345_1;
            double rmu_d2345_2 = r_d2345_2 - mu_d2345_2;
            double rmu_d2345_3 = r_d2345_3 - mu_d2345_3;

            double gmu_d2345_0 = g_d2345_0 - mu_d2345_0;
            double gmu_d2345_1 = g_d2345_1 - mu_d2345_1;
            double gmu_d2345_2 = g_d2345_2 - mu_d2345_2;
            double gmu_d2345_3 = g_d2345_3 - mu_d2345_3;

            double bmu_d2345_0 = b_d2345_0 - mu_d2345_0;
            double bmu_d2345_1 = b_d2345_1 - mu_d2345_1;
            double bmu_d2345_2 = b_d2345_2 - mu_d2345_2;
            double bmu_d2345_3 = b_d2345_3 - mu_d2345_3;


            double rmu2_c2345_0 = rmu_c2345_0 * rmu_c2345_0;
            double rmu2_c2345_1 = rmu_c2345_1 * rmu_c2345_1;
            double rmu2_c2345_2 = rmu_c2345_2 * rmu_c2345_2;
            double rmu2_c2345_3 = rmu_c2345_3 * rmu_c2345_3;

            double gmu2_c2345_0 = gmu_c2345_0 * gmu_c2345_0;
            double gmu2_c2345_1 = gmu_c2345_1 * gmu_c2345_1;
            double gmu2_c2345_2 = gmu_c2345_2 * gmu_c2345_2;
            double gmu2_c2345_3 = gmu_c2345_3 * gmu_c2345_3;

            double bmu2_c2345_0 = bmu_c2345_0 * bmu_c2345_0;
            double bmu2_c2345_1 = bmu_c2345_1 * bmu_c2345_1;
            double bmu2_c2345_2 = bmu_c2345_2 * bmu_c2345_2;
            double bmu2_c2345_3 = bmu_c2345_3 * bmu_c2345_3;

            double rmu2_d2345_0 = rmu_d2345_0 * rmu_d2345_0;
            double rmu2_d2345_1 = rmu_d2345_1 * rmu_d2345_1;
            double rmu2_d2345_2 = rmu_d2345_2 * rmu_d2345_2;
            double rmu2_d2345_3 = rmu_d2345_3 * rmu_d2345_3;

            double gmu2_d2345_0 = gmu_d2345_0 * gmu_d2345_0;
            double gmu2_d2345_1 = gmu_d2345_1 * gmu_d2345_1;
            double gmu2_d2345_2 = gmu_d2345_2 * gmu_d2345_2;
            double gmu2_d2345_3 = gmu_d2345_3 * gmu_d2345_3;

            double bmu2_d2345_0 = bmu_d2345_0 * bmu_d2345_0;
            double bmu2_d2345_1 = bmu_d2345_1 * bmu_d2345_1;
            double bmu2_d2345_2 = bmu_d2345_2 * bmu_d2345_2;
            double bmu2_d2345_3 = bmu_d2345_3 * bmu_d2345_3;


            //sum up rmu2,gmu2,bmu2
            double mu2sum_c2345_0 = rmu2_c2345_0 + gmu2_c2345_0 + bmu2_c2345_0;
            double mu2sum_c2345_1 = rmu2_c2345_1 + gmu2_c2345_1 + bmu2_c2345_1;
            double mu2sum_c2345_2 = rmu2_c2345_2 + gmu2_c2345_2 + bmu2_c2345_2;
            double mu2sum_c2345_3 = rmu2_c2345_3 + gmu2_c2345_3 + bmu2_c2345_3;

            double mu2sum_d2345_0 = rmu2_d2345_0 + gmu2_d2345_0 + bmu2_d2345_0;
            double mu2sum_d2345_1 = rmu2_d2345_1 + gmu2_d2345_1 + bmu2_d2345_1;
            double mu2sum_d2345_2 = rmu2_d2345_2 + gmu2_d2345_2 + bmu2_d2345_2;
            double mu2sum_d2345_3 = rmu2_d2345_3 + gmu2_d2345_3 + bmu2_d2345_3;

            double mu2sumdiv_c2345_0 = mu2sum_c2345_0 / 3.0;
            double mu2sumdiv_c2345_1 = mu2sum_c2345_1 / 3.0;
            double mu2sumdiv_c2345_2 = mu2sum_c2345_2 / 3.0;
            double mu2sumdiv_c2345_3 = mu2sum_c2345_3 / 3.0;

            double mu2sumdiv_d2345_0 = mu2sum_d2345_0 / 3.0;
            double mu2sumdiv_d2345_1 = mu2sum_d2345_1 / 3.0;
            double mu2sumdiv_d2345_2 = mu2sum_d2345_2 / 3.0;
            double mu2sumdiv_d2345_3 = mu2sum_d2345_3 / 3.0;

            double x7_c2345_0 = sqrt(mu2sumdiv_c2345_0);
            double x7_c2345_1 = sqrt(mu2sumdiv_c2345_1);
            double x7_c2345_2 = sqrt(mu2sumdiv_c2345_2);
            double x7_c2345_3 = sqrt(mu2sumdiv_c2345_3);

            double x7_d2345_0 = sqrt(mu2sumdiv_d2345_0);
            double x7_d2345_1 = sqrt(mu2sumdiv_d2345_1);
            double x7_d2345_2 = sqrt(mu2sumdiv_d2345_2);
            double x7_d2345_3 = sqrt(mu2sumdiv_d2345_3);

            double x8_c2345_0 = fabs(x7_c2345_0);
            double x8_c2345_1 = fabs(x7_c2345_1);
            double x8_c2345_2 = fabs(x7_c2345_2);
            double x8_c2345_3 = fabs(x7_c2345_3);

            double x8_d2345_0 = fabs(x7_d2345_0);
            double x8_d2345_1 = fabs(x7_d2345_1);
            double x8_d2345_2 = fabs(x7_d2345_2);
            double x8_d2345_3 = fabs(x7_d2345_3);

            //saturation done, result in x8

            double rz_c2345_0 = r_c2345_0 - 0.5;
            double rz_c2345_1 = r_c2345_1 - 0.5;
            double rz_c2345_2 = r_c2345_2 - 0.5;
            double rz_c2345_3 = r_c2345_3 - 0.5;

            double gz_c2345_0 = g_c2345_0 - 0.5;
            double gz_c2345_1 = g_c2345_1 - 0.5;
            double gz_c2345_2 = g_c2345_2 - 0.5;
            double gz_c2345_3 = g_c2345_3 - 0.5;

            double bz_c2345_0 = b_c2345_0 - 0.5;
            double bz_c2345_1 = b_c2345_1 - 0.5;
            double bz_c2345_2 = b_c2345_2 - 0.5;
            double bz_c2345_3 = b_c2345_3 - 0.5;

            double rz_d2345_0 = r_d2345_0 - 0.5;
            double rz_d2345_1 = r_d2345_1 - 0.5;
            double rz_d2345_2 = r_d2345_2 - 0.5;
            double rz_d2345_3 = r_d2345_3 - 0.5;

            double gz_d2345_0 = g_d2345_0 - 0.5;
            double gz_d2345_1 = g_d2345_1 - 0.5;
            double gz_d2345_2 = g_d2345_2 - 0.5;
            double gz_d2345_3 = g_d2345_3 - 0.5;

            double bz_d2345_0 = b_d2345_0 - 0.5;
            double bz_d2345_1 = b_d2345_1 - 0.5;
            double bz_d2345_2 = b_d2345_2 - 0.5;
            double bz_d2345_3 = b_d2345_3 - 0.5;


            double rz2_c2345_0 = rz_c2345_0 * rz_c2345_0;
            double rz2_c2345_1 = rz_c2345_1 * rz_c2345_1;
            double rz2_c2345_2 = rz_c2345_2 * rz_c2345_2;
            double rz2_c2345_3 = rz_c2345_3 * rz_c2345_3;

            double gz2_c2345_0 = gz_c2345_0 * gz_c2345_0;
            double gz2_c2345_1 = gz_c2345_1 * gz_c2345_1;
            double gz2_c2345_2 = gz_c2345_2 * gz_c2345_2;
            double gz2_c2345_3 = gz_c2345_3 * gz_c2345_3;

            double bz2_c2345_0 = bz_c2345_0 * bz_c2345_0;
            double bz2_c2345_1 = bz_c2345_1 * bz_c2345_1;
            double bz2_c2345_2 = bz_c2345_2 * bz_c2345_2;
            double bz2_c2345_3 = bz_c2345_3 * bz_c2345_3;

            double rz2_d2345_0 = rz_d2345_0 * rz_d2345_0;
            double rz2_d2345_1 = rz_d2345_1 * rz_d2345_1;
            double rz2_d2345_2 = rz_d2345_2 * rz_d2345_2;
            double rz2_d2345_3 = rz_d2345_3 * rz_d2345_3;

            double gz2_d2345_0 = gz_d2345_0 * gz_d2345_0;
            double gz2_d2345_1 = gz_d2345_1 * gz_d2345_1;
            double gz2_d2345_2 = gz_d2345_2 * gz_d2345_2;
            double gz2_d2345_3 = gz_d2345_3 * gz_d2345_3;

            double bz2_d2345_0 = bz_d2345_0 * bz_d2345_0;
            double bz2_d2345_1 = bz_d2345_1 * bz_d2345_1;
            double bz2_d2345_2 = bz_d2345_2 * bz_d2345_2;
            double bz2_d2345_3 = bz_d2345_3 * bz_d2345_3;


            double rexp_c2345_0 = -12.5 * rz2_c2345_0;
            double rexp_c2345_1 = -12.5 * rz2_c2345_1;
            double rexp_c2345_2 = -12.5 * rz2_c2345_2;
            double rexp_c2345_3 = -12.5 * rz2_c2345_3;

            double gexp_c2345_0 = -12.5 * gz2_c2345_0;
            double gexp_c2345_1 = -12.5 * gz2_c2345_1;
            double gexp_c2345_2 = -12.5 * gz2_c2345_2;
            double gexp_c2345_3 = -12.5 * gz2_c2345_3;

            double bexp_c2345_0 = -12.5 * bz2_c2345_0;
            double bexp_c2345_1 = -12.5 * bz2_c2345_1;
            double bexp_c2345_2 = -12.5 * bz2_c2345_2;
            double bexp_c2345_3 = -12.5 * bz2_c2345_3;

            double rexp_d2345_0 = -12.5 * rz2_d2345_0;
            double rexp_d2345_1 = -12.5 * rz2_d2345_1;
            double rexp_d2345_2 = -12.5 * rz2_d2345_2;
            double rexp_d2345_3 = -12.5 * rz2_d2345_3;

            double gexp_d2345_0 = -12.5 * gz2_d2345_0;
            double gexp_d2345_1 = -12.5 * gz2_d2345_1;
            double gexp_d2345_2 = -12.5 * gz2_d2345_2;
            double gexp_d2345_3 = -12.5 * gz2_d2345_3;

            double bexp_d2345_0 = -12.5 * bz2_d2345_0;
            double bexp_d2345_1 = -12.5 * bz2_d2345_1;
            double bexp_d2345_2 = -12.5 * bz2_d2345_2;
            double bexp_d2345_3 = -12.5 * bz2_d2345_3;


            double rexpf_c2345_0 = exp(rexp_c2345_0);
            double rexpf_c2345_1 = exp(rexp_c2345_1);
            double rexpf_c2345_2 = exp(rexp_c2345_2);
            double rexpf_c2345_3 = exp(rexp_c2345_3);

            double gexpf_c2345_0 = exp(gexp_c2345_0);
            double gexpf_c2345_1 = exp(gexp_c2345_1);
            double gexpf_c2345_2 = exp(gexp_c2345_2);
            double gexpf_c2345_3 = exp(gexp_c2345_3);

            double bexpf_c2345_0 = exp(bexp_c2345_0);
            double bexpf_c2345_1 = exp(bexp_c2345_1);
            double bexpf_c2345_2 = exp(bexp_c2345_2);
            double bexpf_c2345_3 = exp(bexp_c2345_3);

            double rexpf_d2345_0 = exp(rexp_d2345_0);
            double rexpf_d2345_1 = exp(rexp_d2345_1);
            double rexpf_d2345_2 = exp(rexp_d2345_2);
            double rexpf_d2345_3 = exp(rexp_d2345_3);

            double gexpf_d2345_0 = exp(gexp_d2345_0);
            double gexpf_d2345_1 = exp(gexp_d2345_1);
            double gexpf_d2345_2 = exp(gexp_d2345_2);
            double gexpf_d2345_3 = exp(gexp_d2345_3);

            double bexpf_d2345_0 = exp(bexp_d2345_0);
            double bexpf_d2345_1 = exp(bexp_d2345_1);
            double bexpf_d2345_2 = exp(bexp_d2345_2);
            double bexpf_d2345_3 = exp(bexp_d2345_3);

            double x9_c2345_0 = rexpf_c2345_0 * gexpf_c2345_0 * bexpf_c2345_0;
            double x9_c2345_1 = rexpf_c2345_1 * gexpf_c2345_1 * bexpf_c2345_1;
            double x9_c2345_2 = rexpf_c2345_2 * gexpf_c2345_2 * bexpf_c2345_2;
            double x9_c2345_3 = rexpf_c2345_3 * gexpf_c2345_3 * bexpf_c2345_3;

            double x9_d2345_0 = rexpf_d2345_0 * gexpf_d2345_0 * bexpf_d2345_0;
            double x9_d2345_1 = rexpf_d2345_1 * gexpf_d2345_1 * bexpf_d2345_1;
            double x9_d2345_2 = rexpf_d2345_2 * gexpf_d2345_2 * bexpf_d2345_2;
            double x9_d2345_3 = rexpf_d2345_3 * gexpf_d2345_3 * bexpf_d2345_3;

            double x10_c2345_0 = fabs(x9_c2345_0);
            double x10_c2345_1 = fabs(x9_c2345_1);
            double x10_c2345_2 = fabs(x9_c2345_2);
            double x10_c2345_3 = fabs(x9_c2345_3);

            double x10_d2345_0 = fabs(x9_d2345_0);
            double x10_d2345_1 = fabs(x9_d2345_1);
            double x10_d2345_2 = fabs(x9_d2345_2);
            double x10_d2345_3 = fabs(x9_d2345_3);

            //wellexposedness done, result in x10

            double x11_c2345_0 = x8_c2345_0 * x10_c2345_0;
            double x11_c2345_1 = x8_c2345_1 * x10_c2345_1;
            double x11_c2345_2 = x8_c2345_2 * x10_c2345_2;
            double x11_c2345_3 = x8_c2345_3 * x10_c2345_3;

            double x11_d2345_0 = x8_d2345_0 * x10_d2345_0;
            double x11_d2345_1 = x8_d2345_1 * x10_d2345_1;
            double x11_d2345_2 = x8_d2345_2 * x10_d2345_2;
            double x11_d2345_3 = x8_d2345_3 * x10_d2345_3;

            //sat,wexp results in x11

            double x12_c2345_0 = x11_c2345_0 * x6_c2345_0;
            double x12_c2345_1 = x11_c2345_1 * x6_c2345_1;
            double x12_c2345_2 = x11_c2345_2 * x6_c2345_2;
            double x12_c2345_3 = x11_c2345_3 * x6_c2345_3;

            double x12_d2345_0 = x11_d2345_0 * x6_d2345_0;
            double x12_d2345_1 = x11_d2345_1 * x6_d2345_1;
            double x12_d2345_2 = x11_d2345_2 * x6_d2345_2;
            double x12_d2345_3 = x11_d2345_3 * x6_d2345_3;

            //sat,wexp,contrast results in x12

            double x13_c2345_0 = x12_c2345_0 + 1.0E-12;
            double x13_c2345_1 = x12_c2345_1 + 1.0E-12;
            double x13_c2345_2 = x12_c2345_2 + 1.0E-12;
            double x13_c2345_3 = x12_c2345_3 + 1.0E-12;

            double x13_d2345_0 = x12_d2345_0 + 1.0E-12;
            double x13_d2345_1 = x12_d2345_1 + 1.0E-12;
            double x13_d2345_2 = x12_d2345_2 + 1.0E-12;
            double x13_d2345_3 = x12_d2345_3 + 1.0E-12;

            dst[c2c] = x13_c2345_0;
            dst[c2c+1] = x13_c2345_1;
            dst[c2c+2] = x13_c2345_2;
            dst[c2c+3] = x13_c2345_3;

            dst[d2c] = x13_d2345_0;
            dst[d2c+1] = x13_d2345_1;
            dst[d2c+2] = x13_d2345_2;
            dst[d2c+3] = x13_d2345_3;

            COST_INC_MUL(160);
            COST_INC_ADD(152);
            COST_INC_DIV(16);
            COST_INC_EXP(24);
            COST_INC_SQRT(8);
            COST_INC_ABS(24);
            COST_INC_LOAD(48);
            COST_INC_STORE(8);
        }
        COST_INC_MUL(12);
        COST_INC_ADD(8);
        COST_INC_LOAD(12);
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
                COST_INC_LOAD(1);
            }
            for (int n = 0; n < nimages; n++){
                int at = (i*c+j);
                double *Wn = W[n];
                Wn[at] = Wn[at] / sum; COST_INC_DIV(1); //beware of division by zero
                COST_INC_LOAD(1);
                COST_INC_STORE(1);
            }
        }
    }
}

#endif
