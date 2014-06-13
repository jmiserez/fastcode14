#ifndef WEIGHTS_C
#define WEIGHTS_C

#include "weights.h"
#include "avx_mathfun.h"

#define USE_AVX_EXP

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

//    double rW = 0.2989;
//    double gW = 0.5870;
//    double bW = 0.1140;

    __m256d rgb0W = _mm256_set_pd (0, 0.1140, 0.5870, 0.2989); //rgb0
    __m256d minus4_256 = _mm256_set1_pd(-4.0);
    __m256d three = _mm256_set1_pd(3.0);
    __m256d onehalf = _mm256_set1_pd(0.5);
    __m256d eps = _mm256_set1_pd(1.0E-12);
    __m256d minustwelvehalf = _mm256_set1_pd(-12.5);
    __m256d signbit_mask256d = _mm256_set1_pd(-0.); // -0. = 1 << 63

    __attribute__ ((aligned (32)))
    __attribute__ ((aligned (32))) double target_c[4];
    __attribute__ ((aligned (32))) double target_d[4];

#ifndef USE_AVX_EXP
    __attribute__ ((aligned (32)))
    __attribute__ ((aligned (32))) double target_rexp_c[4];
    __attribute__ ((aligned (32))) double target_gexp_c[4];
    __attribute__ ((aligned (32))) double target_bexp_c[4];
    __attribute__ ((aligned (32))) double target_rexp_d[4];
    __attribute__ ((aligned (32))) double target_gexp_d[4];
    __attribute__ ((aligned (32))) double target_bexp_d[4];
#endif

    int i = i_start;
    for(i = i_start; i < i_end; i+=2){
        int j = j_start;

        //preload c3456, d3456 for first iteration
        int c1_x3 = 3*(i*c+j-1);
        int d1_x3 = 3*((i+1)*c+j-1);

        //load c3456
//        double c3456_0 = 0;
//        double c3456_1 = 0;
//        double c3456_2 = 0;
//        double c3456_3 = 0;
//        double c3456_4 = 0;
//        double c3456_5 = 0;
//        double c3456_6 = im[c1_x3];
//        double c3456_7 = im[c1_x3+1];
//        double c3456_8 = im[c1_x3+2];
//        double c3456_9 = im[c1_x3+3];
//        double c3456_10 = im[c1_x3+4];
//        double c3456_11 = im[c1_x3+5];

        int64_t leftbitset = ((int64_t)1) << 63;
        __m256i mask0011 = _mm256_setr_epi64x(leftbitset,leftbitset, 0, 0);
        __m256d c3456_4 = _mm256_maskload_pd(&(im[c1_x3-2]),mask0011);
        __m256d c3456_8 = _mm256_loadu_pd(&(im[c1_x3+2]));

        __m128d c3456_6_lo = _mm256_extractf128_pd (c3456_4, 1);// hi
        __m128d c3456_6_hi = _mm256_extractf128_pd (c3456_8, 0);// lo
        __m256d c3456_6 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (c3456_6_lo), c3456_6_hi, 1); //hi

        __m128d c3456_9_tmp = _mm256_extractf128_pd (c3456_8, 1);// hi
        __m128d c3456_9_lo = _mm_shuffle_pd (c3456_6_hi, c3456_9_tmp, 1); //10
        __m128d c3456_9_hi = _mm_permute_pd (c3456_9_tmp, 1); //10
        __m256d c3456_9 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (c3456_9_lo), c3456_9_hi, 1); //hi

        //load d3456
//        double d3456_0 = 0;
//        double d3456_1 = 0;
//        double d3456_2 = 0;
//        double d3456_3 = 0;
//        double d3456_4 = 0;
//        double d3456_5 = 0;
//        double d3456_6 = im[d1_x3];
//        double d3456_7 = im[d1_x3+1];
//        double d3456_8 = im[d1_x3+2];
//        double d3456_9 = im[d1_x3+3];
//        double d3456_10 = im[d1_x3+4];
//        double d3456_11 = im[d1_x3+5];

//        __m25d d3456_0 = _mm256_setzero_pd();
        __m256d d3456_4 = _mm256_maskload_pd(&(im[d1_x3-2]),mask0011);
        __m256d d3456_8 = _mm256_loadu_pd(&(im[d1_x3+2]));

        __m128d d3456_6_lo = _mm256_extractf128_pd (d3456_4, 1);// hi
        __m128d d3456_6_hi = _mm256_extractf128_pd (d3456_8, 0);// lo
        __m256d d3456_6 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (d3456_6_lo), d3456_6_hi, 1); //hi

        __m128d d3456_9_tmp = _mm256_extractf128_pd (d3456_8, 1);// hi
        __m128d d3456_9_lo = _mm_shuffle_pd (d3456_6_hi, d3456_9_tmp, 1); //10
        __m128d d3456_9_hi = _mm_permute_pd (d3456_9_tmp, 1); //10
        __m256d d3456_9 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (d3456_9_lo), d3456_9_hi, 1); //hi

//        double c3456_grey_0 = 0;
//        double c3456_grey_1 = 0;
//        double c3456_grey_2 =
//                rW * c3456_6 +
//                gW * c3456_7 +
//                bW * c3456_8;
//        double c3456_grey_3 =
//                rW * c3456_9 +
//                gW * c3456_10 +
//                bW * c3456_11;

        __m256d zeros = _mm256_setzero_pd();

        __m256d c3456_grey_2 = _mm256_mul_pd(c3456_6,rgb0W);
        __m256d c3456_grey_3 = _mm256_mul_pd(c3456_9,rgb0W);
        __m256d c3456_grey_tmp0 = _mm256_hadd_pd(c3456_grey_2,c3456_grey_3);
        __m256d c3456_grey = _mm256_hadd_pd(zeros,c3456_grey_tmp0);

//        double d3456_grey_0 = 0;
//        double d3456_grey_1 = 0;
//        double d3456_grey_2 =
//                rW * d3456_6 +
//                gW * d3456_7 +
//                bW * d3456_8;
//        double d3456_grey_3 =
//                rW * d3456_9 +
//                gW * d3456_10 +
//                bW * d3456_11;

        __m256d d3456_grey_2 = _mm256_mul_pd(d3456_6,rgb0W);
        __m256d d3456_grey_3 = _mm256_mul_pd(d3456_9,rgb0W);
        __m256d d3456_grey_tmp0 = _mm256_hadd_pd(d3456_grey_2,d3456_grey_3);
        __m256d d3456_grey = _mm256_hadd_pd(zeros,d3456_grey_tmp0);

        for(j = j_start; j < j_end; j+=4){

            //start indices
            int c2c = i*c+j;
            int d2c = (i+1)*c+j;
            int t2_x3 = 3*((i-1)*c+j);
            int c3_x3 = 3*(i*c+j+1);
            int d3_x3 = 3*((i+1)*c+j+1);
            int b2_x3 = 3*((i+2)*c+j);

            //load t2345
//            double t2345_0 = im[t2_x3];
//            double t2345_1 = im[t2_x3+1];
//            double t2345_2 = im[t2_x3+2];
//            double t2345_3 = im[t2_x3+3];
//            double t2345_4 = im[t2_x3+4];
//            double t2345_5 = im[t2_x3+5];
//            double t2345_6 = im[t2_x3+6];
//            double t2345_7 = im[t2_x3+7];
//            double t2345_8 = im[t2_x3+8];
//            double t2345_9 = im[t2_x3+9];
//            double t2345_10 = im[t2_x3+10];
//            double t2345_11 = im[t2_x3+11];

            __m256d t2345_0 = _mm256_loadu_pd(&(im[t2_x3]));
            __m256d t2345_4 = _mm256_loadu_pd(&(im[t2_x3+4]));
            __m256d t2345_8 = _mm256_loadu_pd(&(im[t2_x3+8]));

            __m128d t2345_3_tmp1 = _mm256_extractf128_pd (t2345_4, 0);// lo
            __m128d t2345_3_tmp2 = _mm256_extractf128_pd (t2345_0, 1);// hi
            __m128d t2345_3_lo = _mm_shuffle_pd (t2345_3_tmp1, t2345_3_tmp2, 1); //10
            __m128d t2345_3_hi = _mm_permute_pd (t2345_3_tmp2, 1); //10
            __m256d t2345_3 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (t2345_3_lo), t2345_3_hi, 1); //hi

            __m128d t2345_6_lo = _mm256_extractf128_pd (t2345_4, 1);// hi
            __m128d t2345_6_hi = _mm256_extractf128_pd (t2345_8, 0);// lo
            __m256d t2345_6 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (t2345_6_lo), t2345_6_hi, 1); //hi

            __m128d t2345_9_tmp = _mm256_extractf128_pd (t2345_8, 1);// hi
            __m128d t2345_9_lo = _mm_shuffle_pd (t2345_6_hi, t2345_9_tmp, 1); //10
            __m128d t2345_9_hi = _mm_permute_pd (t2345_9_tmp, 1); //10
            __m256d t2345_9 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (t2345_9_lo), t2345_9_hi, 1); //hi

            //load c12
////            double c12_0 = c3456_6;
////            double c12_1 = c3456_7;
////            double c12_2 = c3456_8;
//            double c12_3 = c3456_9;
//            double c12_4 = c3456_10;
//            double c12_5 = c3456_11;

            __m256d c12_2 = c3456_8;

            //load d12
////            double d12_0 = d3456_6;
////            double d12_1 = d3456_7;
////            double d12_2 = d3456_8;
//            double d12_3 = d3456_9;
//            double d12_4 = d3456_10;
//            double d12_5 = d3456_11;

            __m256d d12_2 = d3456_8;

            //load c3456
//            c3456_0 = im[c3_x3];
//            c3456_1 = im[c3_x3+1];
//            c3456_2 = im[c3_x3+2];
//            c3456_3 = im[c3_x3+3];
//            c3456_4 = im[c3_x3+4];
//            c3456_5 = im[c3_x3+5];
//            c3456_6 = im[c3_x3+6];
//            c3456_7 = im[c3_x3+7];
//            c3456_8 = im[c3_x3+8];
//            c3456_9 = im[c3_x3+9];
//            c3456_10 = im[c3_x3+10];
//            c3456_11 = im[c3_x3+11];

            __m256d c3456_0 = _mm256_loadu_pd(&(im[c3_x3]));
            __m256d c3456_4 = _mm256_loadu_pd(&(im[c3_x3+4]));
            __m256d c3456_8 = _mm256_loadu_pd(&(im[c3_x3+8]));

            __m128d c3456_3_tmp1 = _mm256_extractf128_pd (c3456_4, 0);// lo
            __m128d c3456_3_tmp2 = _mm256_extractf128_pd (c3456_0, 1);// hi
            __m128d c3456_3_lo = _mm_shuffle_pd (c3456_3_tmp1, c3456_3_tmp2, 1); //10
            __m128d c3456_3_hi = _mm_permute_pd (c3456_3_tmp2, 1); //10
            __m256d c3456_3 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (c3456_3_lo), c3456_3_hi, 1); //hi

            __m128d c3456_6_lo = _mm256_extractf128_pd (c3456_4, 1);// hi
            __m128d c3456_6_hi = _mm256_extractf128_pd (c3456_8, 0);// lo
            __m256d c3456_6 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (c3456_6_lo), c3456_6_hi, 1); //hi

            __m128d c3456_9_tmp = _mm256_extractf128_pd (c3456_8, 1);// hi
            __m128d c3456_9_lo = _mm_shuffle_pd (c3456_6_hi, c3456_9_tmp, 1); //10
            __m128d c3456_9_hi = _mm_permute_pd (c3456_9_tmp, 1); //10
            __m256d c3456_9 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (c3456_9_lo), c3456_9_hi, 1); //hi

            //load d3456
//            d3456_0 = im[d3_x3];
//            d3456_1 = im[d3_x3+1];
//            d3456_2 = im[d3_x3+2];
//            d3456_3 = im[d3_x3+3];
//            d3456_4 = im[d3_x3+4];
//            d3456_5 = im[d3_x3+5];
//            d3456_6 = im[d3_x3+6];
//            d3456_7 = im[d3_x3+7];
//            d3456_8 = im[d3_x3+8];
//            d3456_9 = im[d3_x3+9];
//            d3456_10 = im[d3_x3+10];
//            d3456_11 = im[d3_x3+11];

            __m256d d3456_0 = _mm256_loadu_pd(&(im[d3_x3]));
            __m256d d3456_4 = _mm256_loadu_pd(&(im[d3_x3+4]));
            __m256d d3456_8 = _mm256_loadu_pd(&(im[d3_x3+8]));

            __m128d d3456_3_tmp1 = _mm256_extractf128_pd (d3456_4, 0);// lo
            __m128d d3456_3_tmp2 = _mm256_extractf128_pd (d3456_0, 1);// hi
            __m128d d3456_3_lo = _mm_shuffle_pd (d3456_3_tmp1, d3456_3_tmp2, 1); //10
            __m128d d3456_3_hi = _mm_permute_pd (d3456_3_tmp2, 1); //10
            __m256d d3456_3 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (d3456_3_lo), d3456_3_hi, 1); //hi

            __m128d d3456_6_lo = _mm256_extractf128_pd (d3456_4, 1);// hi
            __m128d d3456_6_hi = _mm256_extractf128_pd (d3456_8, 0);// lo
            __m256d d3456_6 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (d3456_6_lo), d3456_6_hi, 1); //hi

            __m128d d3456_9_tmp = _mm256_extractf128_pd (d3456_8, 1);// hi
            __m128d d3456_9_lo = _mm_shuffle_pd (d3456_6_hi, d3456_9_tmp, 1); //10
            __m128d d3456_9_hi = _mm_permute_pd (d3456_9_tmp, 1); //10
            __m256d d3456_9 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (d3456_9_lo), d3456_9_hi, 1); //hi


            //load b2345
//            double b2345_0 = im[b2_x3];
//            double b2345_1 = im[b2_x3+1];
//            double b2345_2 = im[b2_x3+2];
//            double b2345_3 = im[b2_x3+3];
//            double b2345_4 = im[b2_x3+4];
//            double b2345_5 = im[b2_x3+5];
//            double b2345_6 = im[b2_x3+6];
//            double b2345_7 = im[b2_x3+7];
//            double b2345_8 = im[b2_x3+8];
//            double b2345_9 = im[b2_x3+9];
//            double b2345_10 = im[b2_x3+10];
//            double b2345_11 = im[b2_x3+11];

            __m256d b2345_0 = _mm256_loadu_pd(&(im[b2_x3]));
            __m256d b2345_4 = _mm256_loadu_pd(&(im[b2_x3+4]));
            __m256d b2345_8 = _mm256_loadu_pd(&(im[b2_x3+8]));

            __m128d b2345_3_tmp1 = _mm256_extractf128_pd (b2345_4, 0);// lo
            __m128d b2345_3_tmp2 = _mm256_extractf128_pd (b2345_0, 1);// hi
            __m128d b2345_3_lo = _mm_shuffle_pd (b2345_3_tmp1, b2345_3_tmp2, 1); //10
            __m128d b2345_3_hi = _mm_permute_pd (b2345_3_tmp2, 1); //10
            __m256d b2345_3 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (b2345_3_lo), b2345_3_hi, 1); //hi

            __m128d b2345_6_lo = _mm256_extractf128_pd (b2345_4, 1);// hi
            __m128d b2345_6_hi = _mm256_extractf128_pd (b2345_8, 0);// lo
            __m256d b2345_6 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (b2345_6_lo), b2345_6_hi, 1); //hi

            __m128d b2345_9_tmp = _mm256_extractf128_pd (b2345_8, 1);// hi
            __m128d b2345_9_lo = _mm_shuffle_pd (b2345_6_hi, b2345_9_tmp, 1); //10
            __m128d b2345_9_hi = _mm_permute_pd (b2345_9_tmp, 1); //10
            __m256d b2345_9 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (b2345_9_lo), b2345_9_hi, 1); //hi

            //copy old greyscale values
//            double c12_grey_0 = c3456_grey_2;
//            double c12_grey_1 = c3456_grey_3;
//            double d12_grey_0 = d3456_grey_2;
//            double d12_grey_1 = d3456_grey_3;

            __m128d c12_grey = _mm256_extractf128_pd(c3456_grey,1); //extract high pair
            __m128d d12_grey = _mm256_extractf128_pd(d3456_grey,1); //extract high pair

            //calculate greyscales
//            double t2345_grey_0 =
//                    rW * t2345_0 +
//                    gW * t2345_1 +
//                    bW * t2345_2;
//            double t2345_grey_1 =
//                    rW * t2345_3 +
//                    gW * t2345_4 +
//                    bW * t2345_5;
//            double t2345_grey_2 =
//                    rW * t2345_6 +
//                    gW * t2345_7 +
//                    bW * t2345_8;
//            double t2345_grey_3 =
//                    rW * t2345_9 +
//                    gW * t2345_10 +
//                    bW * t2345_11;

            __m256d t2345_grey_0 = _mm256_mul_pd(t2345_0,rgb0W);
            __m256d t2345_grey_1 = _mm256_mul_pd(t2345_3,rgb0W);
            __m256d t2345_grey_2 = _mm256_mul_pd(t2345_6,rgb0W);
            __m256d t2345_grey_3 = _mm256_mul_pd(t2345_9,rgb0W);

            __m256d t2345_grey_tmp1 = _mm256_hadd_pd(t2345_grey_0,t2345_grey_1);
            __m256d t2345_grey_tmp2 = _mm256_hadd_pd(t2345_grey_2,t2345_grey_3);
            __m256d t2345_grey = _mm256_hadd_pd(t2345_grey_tmp1,t2345_grey_tmp2);

//            double b2345_grey_0 =
//                    rW * b2345_0 +
//                    gW * b2345_1 +
//                    bW * b2345_2;
//            double b2345_grey_1 =
//                    rW * b2345_3 +
//                    gW * b2345_4 +
//                    bW * b2345_5;
//            double b2345_grey_2 =
//                    rW * b2345_6 +
//                    gW * b2345_7 +
//                    bW * b2345_8;
//            double b2345_grey_3 =
//                    rW * b2345_9 +
//                    gW * b2345_10 +
//                    bW * b2345_11;

            __m256d b2345_grey_0 = _mm256_mul_pd(b2345_0,rgb0W);
            __m256d b2345_grey_1 = _mm256_mul_pd(b2345_3,rgb0W);
            __m256d b2345_grey_2 = _mm256_mul_pd(b2345_6,rgb0W);
            __m256d b2345_grey_3 = _mm256_mul_pd(b2345_9,rgb0W);

            __m256d b2345_grey_tmp1 = _mm256_hadd_pd(b2345_grey_0,b2345_grey_1);
            __m256d b2345_grey_tmp2 = _mm256_hadd_pd(b2345_grey_2,b2345_grey_3);
            __m256d b2345_grey = _mm256_hadd_pd(b2345_grey_tmp1,b2345_grey_tmp2);

//            c3456_grey_0 =
//                    rW * c3456_0 +
//                    gW * c3456_1 +
//                    bW * c3456_2;
//            c3456_grey_1 =
//                    rW * c3456_3 +
//                    gW * c3456_4 +
//                    bW * c3456_5;
//            c3456_grey_2 =
//                    rW * c3456_6 +
//                    gW * c3456_7 +
//                    bW * c3456_8;
//            c3456_grey_3 =
//                    rW * c3456_9 +
//                    gW * c3456_10 +
//                    bW * c3456_11;

            __m256d c3456_grey_0 = _mm256_mul_pd(c3456_0,rgb0W);
            __m256d c3456_grey_1 = _mm256_mul_pd(c3456_3,rgb0W);
            c3456_grey_2 = _mm256_mul_pd(c3456_6,rgb0W);
            c3456_grey_3 = _mm256_mul_pd(c3456_9,rgb0W);

            __m256d c3456_grey_tmp1 = _mm256_hadd_pd(c3456_grey_0,c3456_grey_1);
            __m256d c3456_grey_tmp2 = _mm256_hadd_pd(c3456_grey_2,c3456_grey_3);
            c3456_grey = _mm256_hadd_pd(c3456_grey_tmp1,c3456_grey_tmp2);

//            d3456_grey_0 =
//                    rW * d3456_0 +
//                    gW * d3456_1 +
//                    bW * d3456_2;
//            d3456_grey_1 =
//                    rW * d3456_3 +
//                    gW * d3456_4 +
//                    bW * d3456_5;
//            d3456_grey_2 =
//                    rW * d3456_6 +
//                    gW * d3456_7 +
//                    bW * d3456_8;
//            d3456_grey_3 =
//                    rW * d3456_9 +
//                    gW * d3456_10 +
//                    bW * d3456_11;

            __m256d d3456_grey_0 = _mm256_mul_pd(d3456_0,rgb0W);
            __m256d d3456_grey_1 = _mm256_mul_pd(d3456_3,rgb0W);
            d3456_grey_2 = _mm256_mul_pd(d3456_6,rgb0W);
            d3456_grey_3 = _mm256_mul_pd(d3456_9,rgb0W);

            __m256d d3456_grey_tmp1 = _mm256_hadd_pd(d3456_grey_0,d3456_grey_1);
            __m256d d3456_grey_tmp2 = _mm256_hadd_pd(d3456_grey_2,d3456_grey_3);
            d3456_grey = _mm256_hadd_pd(d3456_grey_tmp1,d3456_grey_tmp2);

            // shuffle
//            double r_c2345_0 = c12_3;
//            double r_c2345_1 = c3456_0;
//            double r_c2345_2 = c3456_3;
//            double r_c2345_3 = c3456_6;

//            double g_c2345_0 = c12_4;
//            double g_c2345_1 = c3456_1;
//            double g_c2345_2 = c3456_4;
//            double g_c2345_3 = c3456_7;

//            double b_c2345_0 = c12_5;
//            double b_c2345_1 = c3456_2;
//            double b_c2345_2 = c3456_5;
//            double b_c2345_3 = c3456_8;

            __m128d r_c2345_tmp0_1 = _mm256_extractf128_pd (c12_2, 0);// lo
            __m128d r_c2345_tmp0_2 = _mm256_extractf128_pd (c12_2, 1);// hi
            __m128d r_c2345_tmp1_1 = _mm256_extractf128_pd (c3456_0, 0);// lo
            __m128d r_c2345_tmp1_2 = _mm256_extractf128_pd (c3456_0, 1);// hi
            __m128d r_c2345_tmp2_1 = _mm256_extractf128_pd (c3456_3, 0);// lo
            __m128d r_c2345_tmp2_2 = _mm256_extractf128_pd (c3456_3, 1);// hi
            __m128d r_c2345_tmp3_1 = _mm256_extractf128_pd (c3456_6, 0);// lo
            __m128d r_c2345_tmp3_2 = _mm256_extractf128_pd (c3456_6, 1);// hi
            __m128d r_c2345_lo = _mm_shuffle_pd (r_c2345_tmp0_1, r_c2345_tmp1_1, 1); //10
            __m128d r_c2345_hi = _mm_shuffle_pd (r_c2345_tmp2_1, r_c2345_tmp3_1, 0); //00
            __m256d r_c2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (r_c2345_lo), r_c2345_hi, 1); //hi

            __m128d g_c2345_lo = _mm_shuffle_pd (r_c2345_tmp0_2, r_c2345_tmp1_1, 2); //01
            __m128d g_c2345_hi = _mm_shuffle_pd (r_c2345_tmp2_1, r_c2345_tmp3_1, 3); //11
            __m256d g_c2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (g_c2345_lo), g_c2345_hi, 1); //hi

            __m128d b_c2345_lo = _mm_shuffle_pd (r_c2345_tmp0_2, r_c2345_tmp1_2, 1); //10
            __m128d b_c2345_hi = _mm_shuffle_pd (r_c2345_tmp2_2, r_c2345_tmp3_2, 0); //00
            __m256d b_c2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (b_c2345_lo), b_c2345_hi, 1); //hi

//            double r_d2345_0 = d12_3;
//            double r_d2345_1 = d3456_0;
//            double r_d2345_2 = d3456_3;
//            double r_d2345_3 = d3456_6;

//            double g_d2345_0 = d12_4;
//            double g_d2345_1 = d3456_1;
//            double g_d2345_2 = d3456_4;
//            double g_d2345_3 = d3456_7;

//            double b_d2345_0 = d12_5;
//            double b_d2345_1 = d3456_2;
//            double b_d2345_2 = d3456_5;
//            double b_d2345_3 = d3456_8;

            __m128d r_d2345_tmp0_1 = _mm256_extractf128_pd (d12_2, 0);// lo
            __m128d r_d2345_tmp0_2 = _mm256_extractf128_pd (d12_2, 1);// hi
            __m128d r_d2345_tmp1_1 = _mm256_extractf128_pd (d3456_0, 0);// lo
            __m128d r_d2345_tmp1_2 = _mm256_extractf128_pd (d3456_0, 1);// hi
            __m128d r_d2345_tmp2_1 = _mm256_extractf128_pd (d3456_3, 0);// lo
            __m128d r_d2345_tmp2_2 = _mm256_extractf128_pd (d3456_3, 1);// hi
            __m128d r_d2345_tmp3_1 = _mm256_extractf128_pd (d3456_6, 0);// lo
            __m128d r_d2345_tmp3_2 = _mm256_extractf128_pd (d3456_6, 1);// hi
            __m128d r_d2345_lo = _mm_shuffle_pd (r_d2345_tmp0_1, r_d2345_tmp1_1, 1); //10
            __m128d r_d2345_hi = _mm_shuffle_pd (r_d2345_tmp2_1, r_d2345_tmp3_1, 0); //00
            __m256d r_d2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (r_d2345_lo), r_d2345_hi, 1); //hi

            __m128d g_d2345_lo = _mm_shuffle_pd (r_d2345_tmp0_2, r_d2345_tmp1_1, 2); //01
            __m128d g_d2345_hi = _mm_shuffle_pd (r_d2345_tmp2_1, r_d2345_tmp3_1, 3); //11
            __m256d g_d2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (g_d2345_lo), g_d2345_hi, 1); //hi

            __m128d b_d2345_lo = _mm_shuffle_pd (r_d2345_tmp0_2, r_d2345_tmp1_2, 1); //10
            __m128d b_d2345_hi = _mm_shuffle_pd (r_d2345_tmp2_2, r_d2345_tmp3_2, 0); //00
            __m256d b_d2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (b_d2345_lo), b_d2345_hi, 1); //hi


            //shuffle greyscale for outputs
//            double c2345_grey_0 = c12_grey_1;
//            double c2345_grey_1 = c3456_grey_0;
//            double c2345_grey_2 = c3456_grey_1;
//            double c2345_grey_3 = c3456_grey_2;

            __m128d c2345_grey_tmp1 = _mm256_extractf128_pd (c3456_grey, 0);// lo
            __m128d c2345_grey_tmp2 = _mm256_extractf128_pd (c3456_grey, 1);// hi

            __m128d c2345_grey_lo = _mm_shuffle_pd (c12_grey, c2345_grey_tmp1, 1); //10
            __m128d c2345_grey_hi = _mm_shuffle_pd (c2345_grey_tmp1, c2345_grey_tmp2, 1); //10
            __m256d c2345_grey = _mm256_insertf128_pd (_mm256_castpd128_pd256 (c2345_grey_lo), c2345_grey_hi, 1); //hi

//            double d2345_grey_0 = d12_grey_1;
//            double d2345_grey_1 = d3456_grey_0;
//            double d2345_grey_2 = d3456_grey_1;
//            double d2345_grey_3 = d3456_grey_2;

            __m128d d2345_grey_tmp1 = _mm256_extractf128_pd (d3456_grey, 0);// lo
            __m128d d2345_grey_tmp2 = _mm256_extractf128_pd (d3456_grey, 1);// hi

            __m128d d2345_grey_lo = _mm_shuffle_pd (d12_grey, d2345_grey_tmp1, 1); //10
            __m128d d2345_grey_hi = _mm_shuffle_pd (d2345_grey_tmp1, d2345_grey_tmp2, 1); //10
            __m256d d2345_grey = _mm256_insertf128_pd (_mm256_castpd128_pd256 (d2345_grey_lo), d2345_grey_hi, 1); //hi


            //do convolution, center
//            double x1_c2345_0 = c2345_grey_0 * -4.0;
//            double x1_c2345_1 = c2345_grey_1 * -4.0;
//            double x1_c2345_2 = c2345_grey_2 * -4.0;
//            double x1_c2345_3 = c2345_grey_3 * -4.0;

            __m256d x1_c2345 =  _mm256_mul_pd(c2345_grey, minus4_256);

//            double x1_d2345_0 = d2345_grey_0 * -4.0;
//            double x1_d2345_1 = d2345_grey_1 * -4.0;
//            double x1_d2345_2 = d2345_grey_2 * -4.0;
//            double x1_d2345_3 = d2345_grey_3 * -4.0;

            __m256d x1_d2345 =  _mm256_mul_pd(d2345_grey, minus4_256);

            //do convolution, left
//            double x2_c2345_0 = x1_c2345_0 + c12_grey_0;
//            double x2_c2345_1 = x1_c2345_1 + c12_grey_1;
//            double x2_c2345_2 = x1_c2345_2 + c2345_grey_1;
//            double x2_c2345_3 = x1_c2345_3 + c2345_grey_2;

            __m128d x2_c2345_tmp1 = _mm256_extractf128_pd (x1_c2345, 0);// lo
            __m128d x2_c2345_tmp2 = _mm256_extractf128_pd (x1_c2345, 1);// hi
            __m128d x2_c2345_tmp3 = _mm256_extractf128_pd (c2345_grey, 0);// 0
            __m128d x2_c2345_tmp4 = _mm256_extractf128_pd (c2345_grey, 1);// hi

            __m128d x2_c2345_tmp5 = _mm_shuffle_pd (x2_c2345_tmp3, x2_c2345_tmp4, 1); //10
            __m128d x2_c2345_lo = _mm_add_pd(x2_c2345_tmp1,c12_grey);
            __m128d x2_c2345_hi = _mm_add_pd(x2_c2345_tmp2,x2_c2345_tmp5);
            __m256d x2_c2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (x2_c2345_lo), x2_c2345_hi, 1); //hi

//            double x2_d2345_0 = x1_d2345_0 + d12_grey_0;
//            double x2_d2345_1 = x1_d2345_1 + d12_grey_1;
//            double x2_d2345_2 = x1_d2345_2 + d2345_grey_1;
//            double x2_d2345_3 = x1_d2345_3 + d2345_grey_2;

            __m128d x2_d2345_tmp1 = _mm256_extractf128_pd (x1_d2345, 0);// lo
            __m128d x2_d2345_tmp2 = _mm256_extractf128_pd (x1_d2345, 1);// hi
            __m128d x2_d2345_tmp3 = _mm256_extractf128_pd (d2345_grey, 0);// 0
            __m128d x2_d2345_tmp4 = _mm256_extractf128_pd (d2345_grey, 1);// hi

            __m128d x2_d2345_tmp5 = _mm_shuffle_pd (x2_d2345_tmp3, x2_d2345_tmp4, 1); //10
            __m128d x2_d2345_lo = _mm_add_pd(x2_d2345_tmp1,d12_grey);
            __m128d x2_d2345_hi = _mm_add_pd(x2_d2345_tmp2,x2_d2345_tmp5);
            __m256d x2_d2345 = _mm256_insertf128_pd (_mm256_castpd128_pd256 (x2_d2345_lo), x2_d2345_hi, 1); //hi

            //do convolution, right
//            double x3_c2345_0 = x2_c2345_0 + c3456_grey_0;
//            double x3_c2345_1 = x2_c2345_1 + c3456_grey_1;
//            double x3_c2345_2 = x2_c2345_2 + c3456_grey_2;
//            double x3_c2345_3 = x2_c2345_3 + c3456_grey_3;

            __m256d x3_c2345 =  _mm256_add_pd(x2_c2345, c3456_grey);

//            double x3_d2345_0 = x2_d2345_0 + d3456_grey_0;
//            double x3_d2345_1 = x2_d2345_1 + d3456_grey_1;
//            double x3_d2345_2 = x2_d2345_2 + d3456_grey_2;
//            double x3_d2345_3 = x2_d2345_3 + d3456_grey_3;

            __m256d x3_d2345 =  _mm256_add_pd(x2_d2345, d3456_grey);

            //do convolution, top
//            double x4_c2345_0 = x3_c2345_0 + t2345_grey_0;
//            double x4_c2345_1 = x3_c2345_1 + t2345_grey_1;
//            double x4_c2345_2 = x3_c2345_2 + t2345_grey_2;
//            double x4_c2345_3 = x3_c2345_3 + t2345_grey_3;

            __m256d x4_c2345 =  _mm256_add_pd(x3_c2345, t2345_grey);

//            double x4_d2345_0 = x3_d2345_0 + c2345_grey_0;
//            double x4_d2345_1 = x3_d2345_1 + c2345_grey_1;
//            double x4_d2345_2 = x3_d2345_2 + c2345_grey_2;
//            double x4_d2345_3 = x3_d2345_3 + c2345_grey_3;

            __m256d x4_d2345 =  _mm256_add_pd(x3_d2345, c2345_grey);

            //do convolution, bottom
//            double x5_c2345_0 = x4_c2345_0 + d2345_grey_0;
//            double x5_c2345_1 = x4_c2345_1 + d2345_grey_1;
//            double x5_c2345_2 = x4_c2345_2 + d2345_grey_2;
//            double x5_c2345_3 = x4_c2345_3 + d2345_grey_3;

            __m256d x5_c2345 =  _mm256_add_pd(x4_c2345, d2345_grey);

//            double x5_d2345_0 = x4_d2345_0 + b2345_grey_0;
//            double x5_d2345_1 = x4_d2345_1 + b2345_grey_1;
//            double x5_d2345_2 = x4_d2345_2 + b2345_grey_2;
//            double x5_d2345_3 = x4_d2345_3 + b2345_grey_3;

            __m256d x5_d2345 =  _mm256_add_pd(x4_d2345, b2345_grey);

//            double x6_c2345_0 = fabs(x5_c2345_0);
//            double x6_c2345_1 = fabs(x5_c2345_1);
//            double x6_c2345_2 = fabs(x5_c2345_2);
//            double x6_c2345_3 = fabs(x5_c2345_3);

            __m256d x6_c2345 =  _mm256_andnot_pd(x5_c2345, signbit_mask256d);

//            double x6_d2345_0 = fabs(x5_d2345_0);
//            double x6_d2345_1 = fabs(x5_d2345_1);
//            double x6_d2345_2 = fabs(x5_d2345_2);
//            double x6_d2345_3 = fabs(x5_d2345_3);

            __m256d x6_d2345 =  _mm256_andnot_pd(x5_d2345, signbit_mask256d);

            //convolution done, result in x6

//            double rgbsum_c2345_0 = r_c2345_0 + g_c2345_0 + b_c2345_0;
//            double rgbsum_c2345_1 = r_c2345_1 + g_c2345_1 + b_c2345_1;
//            double rgbsum_c2345_2 = r_c2345_2 + g_c2345_2 + b_c2345_2;
//            double rgbsum_c2345_3 = r_c2345_3 + g_c2345_3 + b_c2345_3;

            __m256d rgbsum_c2345_tmp =  _mm256_add_pd(r_c2345, g_c2345);
            __m256d rgbsum_c2345 =  _mm256_add_pd(rgbsum_c2345_tmp, b_c2345);

//            double rgbsum_d2345_0 = r_d2345_0 + g_d2345_0 + b_d2345_0;
//            double rgbsum_d2345_1 = r_d2345_1 + g_d2345_1 + b_d2345_1;
//            double rgbsum_d2345_2 = r_d2345_2 + g_d2345_2 + b_d2345_2;
//            double rgbsum_d2345_3 = r_d2345_3 + g_d2345_3 + b_d2345_3;

            __m256d rgbsum_d2345_tmp =  _mm256_add_pd(r_d2345, g_d2345);
            __m256d rgbsum_d2345 =  _mm256_add_pd(rgbsum_d2345_tmp, b_d2345);

//            double mu_c2345_0 = rgbsum_c2345_0 / 3.0;
//            double mu_c2345_1 = rgbsum_c2345_1 / 3.0;
//            double mu_c2345_2 = rgbsum_c2345_2 / 3.0;
//            double mu_c2345_3 = rgbsum_c2345_3 / 3.0;

            __m256d mu_c2345 =  _mm256_div_pd(rgbsum_c2345, three);

//            double mu_d2345_0 = rgbsum_d2345_0 / 3.0;
//            double mu_d2345_1 = rgbsum_d2345_1 / 3.0;
//            double mu_d2345_2 = rgbsum_d2345_2 / 3.0;
//            double mu_d2345_3 = rgbsum_d2345_3 / 3.0;

            __m256d mu_d2345 =  _mm256_div_pd(rgbsum_d2345, three);

//            double rmu_c2345_0 = r_c2345_0 - mu_c2345_0;
//            double rmu_c2345_1 = r_c2345_1 - mu_c2345_1;
//            double rmu_c2345_2 = r_c2345_2 - mu_c2345_2;
//            double rmu_c2345_3 = r_c2345_3 - mu_c2345_3;

//            double gmu_c2345_0 = g_c2345_0 - mu_c2345_0;
//            double gmu_c2345_1 = g_c2345_1 - mu_c2345_1;
//            double gmu_c2345_2 = g_c2345_2 - mu_c2345_2;
//            double gmu_c2345_3 = g_c2345_3 - mu_c2345_3;

//            double bmu_c2345_0 = b_c2345_0 - mu_c2345_0;
//            double bmu_c2345_1 = b_c2345_1 - mu_c2345_1;
//            double bmu_c2345_2 = b_c2345_2 - mu_c2345_2;
//            double bmu_c2345_3 = b_c2345_3 - mu_c2345_3;

            __m256d rmu_c2345 =  _mm256_sub_pd(r_c2345, mu_c2345);
            __m256d gmu_c2345 =  _mm256_sub_pd(g_c2345, mu_c2345);
            __m256d bmu_c2345 =  _mm256_sub_pd(b_c2345, mu_c2345);


//            double rmu_d2345_0 = r_d2345_0 - mu_d2345_0;
//            double rmu_d2345_1 = r_d2345_1 - mu_d2345_1;
//            double rmu_d2345_2 = r_d2345_2 - mu_d2345_2;
//            double rmu_d2345_3 = r_d2345_3 - mu_d2345_3;

//            double gmu_d2345_0 = g_d2345_0 - mu_d2345_0;
//            double gmu_d2345_1 = g_d2345_1 - mu_d2345_1;
//            double gmu_d2345_2 = g_d2345_2 - mu_d2345_2;
//            double gmu_d2345_3 = g_d2345_3 - mu_d2345_3;

//            double bmu_d2345_0 = b_d2345_0 - mu_d2345_0;
//            double bmu_d2345_1 = b_d2345_1 - mu_d2345_1;
//            double bmu_d2345_2 = b_d2345_2 - mu_d2345_2;
//            double bmu_d2345_3 = b_d2345_3 - mu_d2345_3;

            __m256d rmu_d2345 =  _mm256_sub_pd(r_d2345, mu_d2345);
            __m256d gmu_d2345 =  _mm256_sub_pd(g_d2345, mu_d2345);
            __m256d bmu_d2345 =  _mm256_sub_pd(b_d2345, mu_d2345);


//            double rmu2_c2345_0 = rmu_c2345_0 * rmu_c2345_0;
//            double rmu2_c2345_1 = rmu_c2345_1 * rmu_c2345_1;
//            double rmu2_c2345_2 = rmu_c2345_2 * rmu_c2345_2;
//            double rmu2_c2345_3 = rmu_c2345_3 * rmu_c2345_3;

//            double gmu2_c2345_0 = gmu_c2345_0 * gmu_c2345_0;
//            double gmu2_c2345_1 = gmu_c2345_1 * gmu_c2345_1;
//            double gmu2_c2345_2 = gmu_c2345_2 * gmu_c2345_2;
//            double gmu2_c2345_3 = gmu_c2345_3 * gmu_c2345_3;

//            double bmu2_c2345_0 = bmu_c2345_0 * bmu_c2345_0;
//            double bmu2_c2345_1 = bmu_c2345_1 * bmu_c2345_1;
//            double bmu2_c2345_2 = bmu_c2345_2 * bmu_c2345_2;
//            double bmu2_c2345_3 = bmu_c2345_3 * bmu_c2345_3;

            __m256d rmu2_c2345 =  _mm256_mul_pd(rmu_c2345, rmu_c2345);
            __m256d gmu2_c2345 =  _mm256_mul_pd(gmu_c2345, gmu_c2345);
            __m256d bmu2_c2345 =  _mm256_mul_pd(bmu_c2345, bmu_c2345);


//            double rmu2_d2345_0 = rmu_d2345_0 * rmu_d2345_0;
//            double rmu2_d2345_1 = rmu_d2345_1 * rmu_d2345_1;
//            double rmu2_d2345_2 = rmu_d2345_2 * rmu_d2345_2;
//            double rmu2_d2345_3 = rmu_d2345_3 * rmu_d2345_3;

//            double gmu2_d2345_0 = gmu_d2345_0 * gmu_d2345_0;
//            double gmu2_d2345_1 = gmu_d2345_1 * gmu_d2345_1;
//            double gmu2_d2345_2 = gmu_d2345_2 * gmu_d2345_2;
//            double gmu2_d2345_3 = gmu_d2345_3 * gmu_d2345_3;

//            double bmu2_d2345_0 = bmu_d2345_0 * bmu_d2345_0;
//            double bmu2_d2345_1 = bmu_d2345_1 * bmu_d2345_1;
//            double bmu2_d2345_2 = bmu_d2345_2 * bmu_d2345_2;
//            double bmu2_d2345_3 = bmu_d2345_3 * bmu_d2345_3;

            __m256d rmu2_d2345 =  _mm256_mul_pd(rmu_d2345, rmu_d2345);
            __m256d gmu2_d2345 =  _mm256_mul_pd(gmu_d2345, gmu_d2345);
            __m256d bmu2_d2345 =  _mm256_mul_pd(bmu_d2345, bmu_d2345);

            //sum up rmu2,gmu2,bmu2
//            double mu2sum_c2345_0 = rmu2_c2345_0 + gmu2_c2345_0 + bmu2_c2345_0;
//            double mu2sum_c2345_1 = rmu2_c2345_1 + gmu2_c2345_1 + bmu2_c2345_1;
//            double mu2sum_c2345_2 = rmu2_c2345_2 + gmu2_c2345_2 + bmu2_c2345_2;
//            double mu2sum_c2345_3 = rmu2_c2345_3 + gmu2_c2345_3 + bmu2_c2345_3;

            __m256d mu2sum_c2345_tmp =  _mm256_add_pd(rmu2_c2345, gmu2_c2345);
            __m256d mu2sum_c2345 =  _mm256_add_pd(mu2sum_c2345_tmp, bmu2_c2345);

//            double mu2sum_d2345_0 = rmu2_d2345_0 + gmu2_d2345_0 + bmu2_d2345_0;
//            double mu2sum_d2345_1 = rmu2_d2345_1 + gmu2_d2345_1 + bmu2_d2345_1;
//            double mu2sum_d2345_2 = rmu2_d2345_2 + gmu2_d2345_2 + bmu2_d2345_2;
//            double mu2sum_d2345_3 = rmu2_d2345_3 + gmu2_d2345_3 + bmu2_d2345_3;

            __m256d mu2sum_d2345_tmp =  _mm256_add_pd(rmu2_d2345, gmu2_d2345);
            __m256d mu2sum_d2345 =  _mm256_add_pd(mu2sum_d2345_tmp, bmu2_d2345);

//            double mu2sumdiv_c2345_0 = mu2sum_c2345_0 / 3.0;
//            double mu2sumdiv_c2345_1 = mu2sum_c2345_1 / 3.0;
//            double mu2sumdiv_c2345_2 = mu2sum_c2345_2 / 3.0;
//            double mu2sumdiv_c2345_3 = mu2sum_c2345_3 / 3.0;

            __m256d mu2sumdiv_c2345 =  _mm256_div_pd(mu2sum_c2345, three);

//            double mu2sumdiv_d2345_0 = mu2sum_d2345_0 / 3.0;
//            double mu2sumdiv_d2345_1 = mu2sum_d2345_1 / 3.0;
//            double mu2sumdiv_d2345_2 = mu2sum_d2345_2 / 3.0;
//            double mu2sumdiv_d2345_3 = mu2sum_d2345_3 / 3.0;

            __m256d mu2sumdiv_d2345 =  _mm256_div_pd(mu2sum_d2345, three);

//            double x7_c2345_0 = sqrt(mu2sumdiv_c2345_0);
//            double x7_c2345_1 = sqrt(mu2sumdiv_c2345_1);
//            double x7_c2345_2 = sqrt(mu2sumdiv_c2345_2);
//            double x7_c2345_3 = sqrt(mu2sumdiv_c2345_3);

            __m256d x7_c2345 =  _mm256_sqrt_pd(mu2sumdiv_c2345);

//            double x7_d2345_0 = sqrt(mu2sumdiv_d2345_0);
//            double x7_d2345_1 = sqrt(mu2sumdiv_d2345_1);
//            double x7_d2345_2 = sqrt(mu2sumdiv_d2345_2);
//            double x7_d2345_3 = sqrt(mu2sumdiv_d2345_3);

            __m256d x7_d2345 =  _mm256_sqrt_pd(mu2sumdiv_d2345);

//            double x8_c2345_0 = fabs(x7_c2345_0);
//            double x8_c2345_1 = fabs(x7_c2345_1);
//            double x8_c2345_2 = fabs(x7_c2345_2);
//            double x8_c2345_3 = fabs(x7_c2345_3);

            __m256d x8_c2345 =  _mm256_andnot_pd(x7_c2345, signbit_mask256d);


//            double x8_d2345_0 = fabs(x7_d2345_0);
//            double x8_d2345_1 = fabs(x7_d2345_1);
//            double x8_d2345_2 = fabs(x7_d2345_2);
//            double x8_d2345_3 = fabs(x7_d2345_3);

            __m256d x8_d2345 =  _mm256_andnot_pd(x7_d2345, signbit_mask256d);

            //saturation done, result in x8

//            double rz_c2345_0 = r_c2345_0 - 0.5;
//            double rz_c2345_1 = r_c2345_1 - 0.5;
//            double rz_c2345_2 = r_c2345_2 - 0.5;
//            double rz_c2345_3 = r_c2345_3 - 0.5;

//            double gz_c2345_0 = g_c2345_0 - 0.5;
//            double gz_c2345_1 = g_c2345_1 - 0.5;
//            double gz_c2345_2 = g_c2345_2 - 0.5;
//            double gz_c2345_3 = g_c2345_3 - 0.5;

//            double bz_c2345_0 = b_c2345_0 - 0.5;
//            double bz_c2345_1 = b_c2345_1 - 0.5;
//            double bz_c2345_2 = b_c2345_2 - 0.5;
//            double bz_c2345_3 = b_c2345_3 - 0.5;

            __m256d rz_c2345 =  _mm256_sub_pd(r_c2345, onehalf);
            __m256d gz_c2345 =  _mm256_sub_pd(g_c2345, onehalf);
            __m256d bz_c2345 =  _mm256_sub_pd(b_c2345, onehalf);

//            double rz_d2345_0 = r_d2345_0 - 0.5;
//            double rz_d2345_1 = r_d2345_1 - 0.5;
//            double rz_d2345_2 = r_d2345_2 - 0.5;
//            double rz_d2345_3 = r_d2345_3 - 0.5;

//            double gz_d2345_0 = g_d2345_0 - 0.5;
//            double gz_d2345_1 = g_d2345_1 - 0.5;
//            double gz_d2345_2 = g_d2345_2 - 0.5;
//            double gz_d2345_3 = g_d2345_3 - 0.5;

//            double bz_d2345_0 = b_d2345_0 - 0.5;
//            double bz_d2345_1 = b_d2345_1 - 0.5;
//            double bz_d2345_2 = b_d2345_2 - 0.5;
//            double bz_d2345_3 = b_d2345_3 - 0.5;

            __m256d rz_d2345 =  _mm256_sub_pd(r_d2345, onehalf);
            __m256d gz_d2345 =  _mm256_sub_pd(g_d2345, onehalf);
            __m256d bz_d2345 =  _mm256_sub_pd(b_d2345, onehalf);

//            double rz2_c2345_0 = rz_c2345_0 * rz_c2345_0;
//            double rz2_c2345_1 = rz_c2345_1 * rz_c2345_1;
//            double rz2_c2345_2 = rz_c2345_2 * rz_c2345_2;
//            double rz2_c2345_3 = rz_c2345_3 * rz_c2345_3;

//            double gz2_c2345_0 = gz_c2345_0 * gz_c2345_0;
//            double gz2_c2345_1 = gz_c2345_1 * gz_c2345_1;
//            double gz2_c2345_2 = gz_c2345_2 * gz_c2345_2;
//            double gz2_c2345_3 = gz_c2345_3 * gz_c2345_3;

//            double bz2_c2345_0 = bz_c2345_0 * bz_c2345_0;
//            double bz2_c2345_1 = bz_c2345_1 * bz_c2345_1;
//            double bz2_c2345_2 = bz_c2345_2 * bz_c2345_2;
//            double bz2_c2345_3 = bz_c2345_3 * bz_c2345_3;

            __m256d rz2_c2345 =  _mm256_mul_pd(rz_c2345, rz_c2345);
            __m256d gz2_c2345 =  _mm256_mul_pd(gz_c2345, gz_c2345);
            __m256d bz2_c2345 =  _mm256_mul_pd(bz_c2345, bz_c2345);

//            double rz2_d2345_0 = rz_d2345_0 * rz_d2345_0;
//            double rz2_d2345_1 = rz_d2345_1 * rz_d2345_1;
//            double rz2_d2345_2 = rz_d2345_2 * rz_d2345_2;
//            double rz2_d2345_3 = rz_d2345_3 * rz_d2345_3;

//            double gz2_d2345_0 = gz_d2345_0 * gz_d2345_0;
//            double gz2_d2345_1 = gz_d2345_1 * gz_d2345_1;
//            double gz2_d2345_2 = gz_d2345_2 * gz_d2345_2;
//            double gz2_d2345_3 = gz_d2345_3 * gz_d2345_3;

//            double bz2_d2345_0 = bz_d2345_0 * bz_d2345_0;
//            double bz2_d2345_1 = bz_d2345_1 * bz_d2345_1;
//            double bz2_d2345_2 = bz_d2345_2 * bz_d2345_2;
//            double bz2_d2345_3 = bz_d2345_3 * bz_d2345_3;

            __m256d rz2_d2345 =  _mm256_mul_pd(rz_d2345, rz_d2345);
            __m256d gz2_d2345 =  _mm256_mul_pd(gz_d2345, gz_d2345);
            __m256d bz2_d2345 =  _mm256_mul_pd(bz_d2345, bz_d2345);

//            double rexp_c2345_0 = -12.5 * rz2_c2345_0;
//            double rexp_c2345_1 = -12.5 * rz2_c2345_1;
//            double rexp_c2345_2 = -12.5 * rz2_c2345_2;
//            double rexp_c2345_3 = -12.5 * rz2_c2345_3;

//            double gexp_c2345_0 = -12.5 * gz2_c2345_0;
//            double gexp_c2345_1 = -12.5 * gz2_c2345_1;
//            double gexp_c2345_2 = -12.5 * gz2_c2345_2;
//            double gexp_c2345_3 = -12.5 * gz2_c2345_3;

//            double bexp_c2345_0 = -12.5 * bz2_c2345_0;
//            double bexp_c2345_1 = -12.5 * bz2_c2345_1;
//            double bexp_c2345_2 = -12.5 * bz2_c2345_2;
//            double bexp_c2345_3 = -12.5 * bz2_c2345_3;

            __m256d rexp_c2345 =  _mm256_mul_pd(rz2_c2345, minustwelvehalf);
            __m256d gexp_c2345 =  _mm256_mul_pd(gz2_c2345, minustwelvehalf);
            __m256d bexp_c2345 =  _mm256_mul_pd(bz2_c2345, minustwelvehalf);

//            double rexp_d2345_0 = -12.5 * rz2_d2345_0;
//            double rexp_d2345_1 = -12.5 * rz2_d2345_1;
//            double rexp_d2345_2 = -12.5 * rz2_d2345_2;
//            double rexp_d2345_3 = -12.5 * rz2_d2345_3;

//            double gexp_d2345_0 = -12.5 * gz2_d2345_0;
//            double gexp_d2345_1 = -12.5 * gz2_d2345_1;
//            double gexp_d2345_2 = -12.5 * gz2_d2345_2;
//            double gexp_d2345_3 = -12.5 * gz2_d2345_3;

//            double bexp_d2345_0 = -12.5 * bz2_d2345_0;
//            double bexp_d2345_1 = -12.5 * bz2_d2345_1;
//            double bexp_d2345_2 = -12.5 * bz2_d2345_2;
//            double bexp_d2345_3 = -12.5 * bz2_d2345_3;

            __m256d rexp_d2345 =  _mm256_mul_pd(rz2_d2345, minustwelvehalf);
            __m256d gexp_d2345 =  _mm256_mul_pd(gz2_d2345, minustwelvehalf);
            __m256d bexp_d2345 =  _mm256_mul_pd(bz2_d2345, minustwelvehalf);

//            double rexpf_c2345_0 = exp(rexp_c2345_0);
//            double rexpf_c2345_1 = exp(rexp_c2345_1);
//            double rexpf_c2345_2 = exp(rexp_c2345_2);
//            double rexpf_c2345_3 = exp(rexp_c2345_3);

//            double gexpf_c2345_0 = exp(gexp_c2345_0);
//            double gexpf_c2345_1 = exp(gexp_c2345_1);
//            double gexpf_c2345_2 = exp(gexp_c2345_2);
//            double gexpf_c2345_3 = exp(gexp_c2345_3);

//            double bexpf_c2345_0 = exp(bexp_c2345_0);
//            double bexpf_c2345_1 = exp(bexp_c2345_1);
//            double bexpf_c2345_2 = exp(bexp_c2345_2);
//            double bexpf_c2345_3 = exp(bexp_c2345_3);

//            double rexpf_d2345_0 = exp(rexp_d2345_0);
//            double rexpf_d2345_1 = exp(rexp_d2345_1);
//            double rexpf_d2345_2 = exp(rexp_d2345_2);
//            double rexpf_d2345_3 = exp(rexp_d2345_3);

//            double gexpf_d2345_0 = exp(gexp_d2345_0);
//            double gexpf_d2345_1 = exp(gexp_d2345_1);
//            double gexpf_d2345_2 = exp(gexp_d2345_2);
//            double gexpf_d2345_3 = exp(gexp_d2345_3);

//            double bexpf_d2345_0 = exp(bexp_d2345_0);
//            double bexpf_d2345_1 = exp(bexp_d2345_1);
//            double bexpf_d2345_2 = exp(bexp_d2345_2);
//            double bexpf_d2345_3 = exp(bexp_d2345_3);


#ifdef USE_AVX_EXP
            //using avx_fun.h library here, unfortunately this is single precision only
            //casts do not generate code
            __m256d rexpf_c2345 =  _mm256_castps_pd(exp256_ps(_mm256_castpd_ps (rexp_c2345)));
            __m256d gexpf_c2345 =  _mm256_castps_pd(exp256_ps(_mm256_castpd_ps (gexp_c2345)));
            __m256d bexpf_c2345 =  _mm256_castps_pd(exp256_ps(_mm256_castpd_ps (bexp_c2345)));

            __m256d rexpf_d2345 =  _mm256_castps_pd(exp256_ps(_mm256_castpd_ps (rexp_d2345)));
            __m256d gexpf_d2345 =  _mm256_castps_pd(exp256_ps(_mm256_castpd_ps (gexp_d2345)));
            __m256d bexpf_d2345 =  _mm256_castps_pd(exp256_ps(_mm256_castpd_ps (bexp_d2345)));
#else
            //This is very slow. Most of the performance gain from using AVX disappears if we do this.
            _mm256_store_pd(&(target_rexp_c[0]),rexp_c2345);
            _mm256_store_pd(&(target_gexp_c[0]),gexp_c2345);
            _mm256_store_pd(&(target_bexp_c[0]),bexp_c2345);
            _mm256_store_pd(&(target_rexp_d[0]),rexp_d2345);
            _mm256_store_pd(&(target_gexp_d[0]),gexp_d2345);
            _mm256_store_pd(&(target_bexp_d[0]),bexp_d2345);

            target_rexp_c[0] = exp(target_rexp_c[0]);
            target_rexp_c[1] = exp(target_rexp_c[1]);
            target_rexp_c[2] = exp(target_rexp_c[2]);
            target_rexp_c[3] = exp(target_rexp_c[3]);

            target_gexp_c[0] = exp(target_gexp_c[0]);
            target_gexp_c[1] = exp(target_gexp_c[1]);
            target_gexp_c[2] = exp(target_gexp_c[2]);
            target_gexp_c[3] = exp(target_gexp_c[3]);

            target_bexp_c[0] = exp(target_bexp_c[0]);
            target_bexp_c[1] = exp(target_bexp_c[1]);
            target_bexp_c[2] = exp(target_bexp_c[2]);
            target_bexp_c[3] = exp(target_bexp_c[3]);

            target_rexp_d[0] = exp(target_rexp_d[0]);
            target_rexp_d[1] = exp(target_rexp_d[1]);
            target_rexp_d[2] = exp(target_rexp_d[2]);
            target_rexp_d[3] = exp(target_rexp_d[3]);

            target_gexp_d[0] = exp(target_gexp_d[0]);
            target_gexp_d[1] = exp(target_gexp_d[1]);
            target_gexp_d[2] = exp(target_gexp_d[2]);
            target_gexp_d[3] = exp(target_gexp_d[3]);

            target_bexp_d[0] = exp(target_bexp_d[0]);
            target_bexp_d[1] = exp(target_bexp_d[1]);
            target_bexp_d[2] = exp(target_bexp_d[2]);
            target_bexp_d[3] = exp(target_bexp_d[3]);

            __m256d rexpf_c2345 = _mm256_load_pd(&(target_rexp_c[0]));
            __m256d gexpf_c2345 = _mm256_load_pd(&(target_gexp_c[0]));
            __m256d bexpf_c2345 = _mm256_load_pd(&(target_bexp_c[0]));
            __m256d rexpf_d2345 = _mm256_load_pd(&(target_rexp_d[0]));
            __m256d gexpf_d2345 = _mm256_load_pd(&(target_gexp_d[0]));
            __m256d bexpf_d2345 = _mm256_load_pd(&(target_bexp_d[0]));

#endif

//            double x9_c2345_0 = rexpf_c2345_0 * gexpf_c2345_0 * bexpf_c2345_0;
//            double x9_c2345_1 = rexpf_c2345_1 * gexpf_c2345_1 * bexpf_c2345_1;
//            double x9_c2345_2 = rexpf_c2345_2 * gexpf_c2345_2 * bexpf_c2345_2;
//            double x9_c2345_3 = rexpf_c2345_3 * gexpf_c2345_3 * bexpf_c2345_3;

            __m256d x9_c2345_tmp =  _mm256_mul_pd(rexpf_c2345, gexpf_c2345);
            __m256d x9_c2345 =  _mm256_mul_pd(x9_c2345_tmp, bexpf_c2345);

//            double x9_d2345_0 = rexpf_d2345_0 * gexpf_d2345_0 * bexpf_d2345_0;
//            double x9_d2345_1 = rexpf_d2345_1 * gexpf_d2345_1 * bexpf_d2345_1;
//            double x9_d2345_2 = rexpf_d2345_2 * gexpf_d2345_2 * bexpf_d2345_2;
//            double x9_d2345_3 = rexpf_d2345_3 * gexpf_d2345_3 * bexpf_d2345_3;

            __m256d x9_d2345_tmp =  _mm256_mul_pd(rexpf_d2345, gexpf_d2345);
            __m256d x9_d2345 =  _mm256_mul_pd(x9_d2345_tmp, bexpf_d2345);

//            double x10_c2345_0 = fabs(x9_c2345_0);
//            double x10_c2345_1 = fabs(x9_c2345_1);
//            double x10_c2345_2 = fabs(x9_c2345_2);
//            double x10_c2345_3 = fabs(x9_c2345_3);

            __m256d x10_c2345 =  _mm256_andnot_pd(x9_c2345, signbit_mask256d);

//            double x10_d2345_0 = fabs(x9_d2345_0);
//            double x10_d2345_1 = fabs(x9_d2345_1);
//            double x10_d2345_2 = fabs(x9_d2345_2);
//            double x10_d2345_3 = fabs(x9_d2345_3);

            __m256d x10_d2345 =  _mm256_andnot_pd(x9_d2345, signbit_mask256d);

            //wellexposedness done, result in x10

//            double x11_c2345_0 = x8_c2345_0 * x10_c2345_0;
//            double x11_c2345_1 = x8_c2345_1 * x10_c2345_1;
//            double x11_c2345_2 = x8_c2345_2 * x10_c2345_2;
//            double x11_c2345_3 = x8_c2345_3 * x10_c2345_3;

            __m256d x11_c2345 =  _mm256_mul_pd(x8_c2345, x10_c2345);

//            double x11_d2345_0 = x8_d2345_0 * x10_d2345_0;
//            double x11_d2345_1 = x8_d2345_1 * x10_d2345_1;
//            double x11_d2345_2 = x8_d2345_2 * x10_d2345_2;
//            double x11_d2345_3 = x8_d2345_3 * x10_d2345_3;

            __m256d x11_d2345 =  _mm256_mul_pd(x8_d2345, x10_d2345);

            //sat,wexp results in x11

//            double x12_c2345_0 = x11_c2345_0 * x6_c2345_0;
//            double x12_c2345_1 = x11_c2345_1 * x6_c2345_1;
//            double x12_c2345_2 = x11_c2345_2 * x6_c2345_2;
//            double x12_c2345_3 = x11_c2345_3 * x6_c2345_3;

            __m256d x12_c2345 =  _mm256_mul_pd(x11_c2345, x6_c2345);

//            double x12_d2345_0 = x11_d2345_0 * x6_d2345_0;
//            double x12_d2345_1 = x11_d2345_1 * x6_d2345_1;
//            double x12_d2345_2 = x11_d2345_2 * x6_d2345_2;
//            double x12_d2345_3 = x11_d2345_3 * x6_d2345_3;

            __m256d x12_d2345 =  _mm256_mul_pd(x11_d2345, x6_d2345);

            //sat,wexp,contrast results in x12

//            double x13_c2345_0 = x12_c2345_0 + 1.0E-12;
//            double x13_c2345_1 = x12_c2345_1 + 1.0E-12;
//            double x13_c2345_2 = x12_c2345_2 + 1.0E-12;
//            double x13_c2345_3 = x12_c2345_3 + 1.0E-12;

            __m256d x13_c2345 =  _mm256_add_pd(x12_c2345, eps);

//            double x13_d2345_0 = x12_d2345_0 + 1.0E-12;
//            double x13_d2345_1 = x12_d2345_1 + 1.0E-12;
//            double x13_d2345_2 = x12_d2345_2 + 1.0E-12;
//            double x13_d2345_3 = x12_d2345_3 + 1.0E-12;

            __m256d x13_d2345 =  _mm256_add_pd(x12_d2345, eps);


            //Store 8 results at once
            _mm256_store_pd(&(target_c[0]),x13_c2345);
            _mm256_store_pd(&(target_d[0]),x13_d2345);
            memcpy(&(dst[c2c]),&(target_c[0]),4*sizeof(double));
            memcpy(&(dst[d2c]),&(target_d[0]),4*sizeof(double));

            COST_INC_MUL(160);
            COST_INC_ADD(144);
            COST_INC_DIV(16);
            COST_INC_EXP(24);
            COST_INC_SQRT(8);
            COST_INC_ABS(24);
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