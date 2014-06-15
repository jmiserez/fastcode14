#ifndef WEIGHTS_C
#define WEIGHTS_C

#include "weights.h"
#define USE_SSE2
#include "sse_mathfun.h"

/**
  * Calculate all values in one step per pixel. Requires grabbing the neighboring pixels.
  */
FORCE_INLINE double single_pixel(
        double *im, int center, int top, int left, int right, int bottom,
        const __m256i mask1110,
        const __m256d rgb0W,
        const __m256d onehalf,
        const __m256d minustwelvehalf){
//    double r = im[center];
//    double g = im[center+1];
//    double b = im[center+2];

//    double r1 = im[top];
//    double g1 = im[top+1];
//    double b1 = im[top+2];
//    double r2 = im[left];
//    double g2 = im[left+1];
//    double b2 = im[left+2];
//    double r3 = im[right];
//    double g3 = im[right+1];
//    double b3 = im[right+2];
//    double r4 = im[bottom];
//    double g4 = im[bottom+1];
//    double b4 = im[bottom+2];

    __m256d c = _mm256_maskload_pd(&(im[center]),mask1110);
    __m256d c1 = _mm256_loadu_pd(&(im[top]));
    __m256d c2 = _mm256_loadu_pd(&(im[left]));
    __m256d c3 = _mm256_loadu_pd(&(im[right]));
    __m256d c4 = _mm256_loadu_pd(&(im[bottom]));

    COST_INC_LOAD(20);

//    double grey = rw * r + gw * g + bw * b;
//    double grey1 = rw * r1 + gw * g1 + bw * b1;
//    double grey2 = rw * r2 + gw * g2 + bw * b2;
//    double grey3 = rw * r3 + gw * g3 + bw * b3;
//    double grey4 = rw * r4 + gw * g4 + bw * b4;

    __m256d greyc = _mm256_mul_pd(c,rgb0W);
    __m256d grey1 = _mm256_mul_pd(c1,rgb0W);
    __m256d grey2 = _mm256_mul_pd(c2,rgb0W);
    __m256d grey3 = _mm256_mul_pd(c3,rgb0W);
    __m256d grey4 = _mm256_mul_pd(c4,rgb0W);

    //AVX: double: horizontal add for 1 vector
     __m256d c_perm = _mm256_permute2f128_pd(c, c, 0b00100001);//1,2
     __m256d c_h = _mm256_hadd_pd(c,c_perm);
     __m128d c_h_lo = _mm256_extractf128_pd (c_h, 0);// lo
     __m128d c_h_hi = _mm256_extractf128_pd (c_h, 1);// hi
     double c_hsum_lo = _mm_cvtsd_f64(c_h_lo);
     double c_hsum_hi = _mm_cvtsd_f64(c_h_hi);
     double c_hsum = c_hsum_lo + c_hsum_hi;

     //AVX: double: horizontal add for 1 vector
      __m256d greyc_perm = _mm256_permute2f128_pd(greyc, greyc, 0b00100001);//1,2
      __m256d greyc_h = _mm256_hadd_pd(greyc,greyc_perm);
      __m128d greyc_h_lo = _mm256_extractf128_pd (greyc_h, 0);// lo
      __m128d greyc_h_hi = _mm256_extractf128_pd (greyc_h, 1);// hi
      double greyc_hsum_lo = _mm_cvtsd_f64(greyc_h_lo);
      double greyc_hsum_hi = _mm_cvtsd_f64(greyc_h_hi);
      double greyc_hsum = greyc_hsum_lo + greyc_hsum_hi;

    //AVX: _m256d: horizontal add for 4 vectors at once
    __m256d grey12 = _mm256_hadd_pd(grey1,grey2);
    __m256d grey34 = _mm256_hadd_pd(grey3,grey4);
    __m256d grey_1234_blend = _mm256_blend_pd(grey12, grey34, 0b1100); //0011
    __m256d grey_1234_perm = _mm256_permute2f128_pd(grey12, grey34, 0b00100001);//1,2
    __m256d grey_1234 =  _mm256_add_pd(grey_1234_perm, grey_1234_blend);

    //AVX: double: horizontal add for 1 vector
     __m256d grey1234_perm = _mm256_permute2f128_pd(grey_1234, grey_1234, 0b00100001);//1,2
     __m256d grey1234_h = _mm256_hadd_pd(grey_1234,grey1234_perm);
     __m128d grey1234_h_lo = _mm256_extractf128_pd (grey1234_h, 0);// lo
     __m128d grey1234_h_hi = _mm256_extractf128_pd (grey1234_h, 1);// hi
     double grey1234_hsum_lo = _mm_cvtsd_f64(grey1234_h_lo);
     double grey1234_hsum_hi = _mm_cvtsd_f64(grey1234_h_hi);
     double grey1234_sum = grey1234_hsum_lo + grey1234_hsum_hi;

    COST_INC_ADD(10); //+ operations wasted on AVX
    COST_INC_MUL(15); //+ operations wasted on AVX

    double mu = c_hsum / 3.0;
    COST_INC_ADD(2);
    COST_INC_DIV(1);

//    double rmu = r-mu;
//    double gmu = g-mu;
//    double bmu = b-mu;

    __m256d c_mu = _mm256_set1_pd(mu);
    __m256d c_rgbmu = _mm256_sub_pd(c,c_mu);
    COST_INC_ADD(3); //+1 operations wasted on AVX

//    double rz = r-0.5;
//    double gz = g-0.5;
//    double bz = b-0.5;

    __m256d c_rgbz = _mm256_sub_pd(c,onehalf);
    COST_INC_ADD(3); //+1 operations wasted on AVX

//    double rzrz = rz*rz;
//    double gzgz = gz*gz;
//    double bzbz = bz*bz;

    __m256d c_rgbz_sq = _mm256_mul_pd(c_rgbz,c_rgbz);
    COST_INC_MUL(3); //+1 operations wasted on AVX

//    double re = exp(-12.5*rzrz);
//    double ge = exp(-12.5*gzgz);
//    double be = exp(-12.5*bzbz);

    __m256d c_rgbe_tmp = _mm256_mul_pd(minustwelvehalf,c_rgbz_sq);

    __m128 c_rgbe_tmp_ps = _mm256_cvtpd_ps(c_rgbe_tmp);
    __m128 c_rgbe_ps = exp_ps(c_rgbe_tmp_ps);
    __m256d c_rgbe = _mm256_cvtps_pd(c_rgbe_ps);

    COST_INC_EXP(3);
    COST_INC_MUL(3); //+1 operations wasted on AVX

//    double t1 = sqrt((rmu*rmu + gmu*gmu + bmu*bmu)/3.0);
    __m256d c_rgbmu_sq = _mm256_mul_pd(c_rgbmu,c_rgbmu);

    __m128d t1_tmp1_lo = _mm256_extractf128_pd (c_rgbmu_sq, 0);// lo
    __m128d t1_tmp1_hi = _mm256_extractf128_pd (c_rgbmu_sq, 1);// hi
    __m128d t1_tmp1_lo_sum = _mm_hadd_pd (t1_tmp1_lo, t1_tmp1_lo);
    double t1_tmp1_hi_lo = _mm_cvtsd_f64(t1_tmp1_hi);
    double t1_tmp1_lo_sum_lo = _mm_cvtsd_f64(t1_tmp1_lo_sum);

    double t1_tmp1 = t1_tmp1_lo_sum_lo + t1_tmp1_hi_lo;

    double t1_tmp2 = t1_tmp1 / 3.0;
    double t1 = sqrt(t1_tmp2);

    COST_INC_SQRT(1);
    COST_INC_ADD(3);
    COST_INC_MUL(3); //+1 operations wasted on AVX
    COST_INC_DIV(1);
    double t2 = fabs(t1);
    COST_INC_ABS(1);

//    double t3 = re*ge*be;

    __m128d t3_tmp1_lo = _mm256_extractf128_pd (c_rgbe, 0);// lo
    __m128d t3_tmp1_hi = _mm256_extractf128_pd (c_rgbe, 1);// hi

    double t3_tmp1_lo_lo = _mm_cvtsd_f64(t3_tmp1_lo);
    double t3_tmp1_hi_lo = _mm_cvtsd_f64(t3_tmp1_hi);
    __m128d t3_tmp1_lo_swapped = _mm_permute_pd(t3_tmp1_lo, 1);// swap
    double t3_tmp1_lo_hi = _mm_cvtsd_f64(t3_tmp1_lo_swapped);

    double t3 = t3_tmp1_lo_lo * t3_tmp1_lo_hi * t3_tmp1_hi_lo;

    COST_INC_MUL(2);
    double t4 = fabs(t3);
    COST_INC_ABS(1);

    double t5 = t2 * t4;
    COST_INC_MUL(1);

//    double t6 = -4.0*grey+grey1+grey2+grey3+grey4;

    double minusfour_times_grey = -4.0*greyc_hsum;
    double t6 = minusfour_times_grey+grey1234_sum;

    COST_INC_MUL(1);
    COST_INC_ADD(2); //2 operations saved due to AVX

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
    const double rw = 0.2989;
    const double gw = 0.5870;
    const double bw = 0.1140;
    const __m256d rgb0W = _mm256_setr_pd (rw,gw,bw,0); //rgb0
    const __m256d onehalf = _mm256_set1_pd(0.5);
    const __m256d minustwelvehalf = _mm256_set1_pd(-12.5);
    int64_t hi_bit_set = ((int64_t)1) << 63;
    const __m256i mask1110 = _mm256_setr_epi64x(hi_bit_set,hi_bit_set, hi_bit_set, 0);
    for(int i = 1; i < r-1; i++){
        for(int j = 1; j < c-1; j++){
            int center =  i*c+j;
            int top    = (i-1)*c+j;
            int bottom = (i+1)*c+j;
            int left   = i*c+j-1;
            int right   = i*c+j+1;
            dst[center] = single_pixel(im,
                                       3*center,
                                       3*top,3*left,3*right,3*bottom,
                                       mask1110,rgb0W,onehalf,minustwelvehalf);
            COST_INC_STORE(1);
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
                                   3*center,3*left,3*right,3*bottom,
                                   mask1110,rgb0W,onehalf,minustwelvehalf);
        COST_INC_STORE(1);
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
                                   3*top,3*left,3*right,3*center,
                                   mask1110,rgb0W,onehalf,minustwelvehalf);
        COST_INC_STORE(1);
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
                                   3*top,3*center,3*right,3*bottom,
                                   mask1110,rgb0W,onehalf,minustwelvehalf);
        COST_INC_STORE(1);
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
                                   3*top,3*left,3*center,3*bottom,
                                   mask1110,rgb0W,onehalf,minustwelvehalf);
        COST_INC_STORE(1);
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
                               3*center,3*center,3*right,3*bottom,
                               mask1110,rgb0W,onehalf,minustwelvehalf);
    COST_INC_STORE(1);
    //top right
    i = 0;
    j = c-1;
    center =  i*c+j;
    bottom = (i+1)*c+j;
    int left   = i*c+j-1;
    dst[center] = single_pixel(im,
                               3*center,
                               3*center,3*left,3*center,3*bottom,
                               mask1110,rgb0W,onehalf,minustwelvehalf);
    COST_INC_STORE(1);
    //bottom left
    i = r-1;
    j = 0;
    center =  i*c+j;
    int top    = (i-1)*c+j;
    right   = i*c+j+1;
    dst[center] = single_pixel(im,
                               3*center,
                               3*top,3*center,3*right,3*center,
                               mask1110,rgb0W,onehalf,minustwelvehalf);
    COST_INC_STORE(1);
    //bottom right
    i = r-1;
    j = c-1;
    center =  i*c+j;
    top    = (i-1)*c+j;
    left   = i*c+j-1;
    dst[center] = single_pixel(im,
                               3*center,
                               3*top,3*left,3*center,3*center,
                               mask1110,rgb0W,onehalf,minustwelvehalf);
    COST_INC_STORE(1);
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
