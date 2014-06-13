#include "fusion_includes.h"

FORCE_INLINE void _weights_inner(double *im, double **imC, int at, int real_at, int top, int left, int right, int bottom, double *Wn){
    double r = im[at];
    double g = im[at+1];
    double b = im[at+2];

    //store individual channels
    //TODO SIMD nocache
    imC[0][real_at] = r;
    imC[1][real_at] = g;
    imC[2][real_at] = b;

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
    COST_INC_ADD(2);
    COST_INC_MUL(3);
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

    Wn[real_at] = t9;
}

FORCE_INLINE void _weights_fallback(
        int rSTART, int rEND, int cSTART, int cEND,
        double *im, double **imC, double *Wn,
        uint32_t r, uint32_t c){
    //center rows
    for(int i = rSTART; i < rEND; i++){
        int itop = i-1;
        int ibottom = i+1;
        if(i == 0){
            itop = 0;
        } else if(i == r-1){
            ibottom = r-1;
        }
        //center cols
        for(int j = cSTART; j < cEND; j++){
            int jleft = j-1;
            int jright = j+1;
            if(j == 0){
                jleft = 0;
            } else if(j == c-1){
                jright = c-1;
            }
            int real_at = (i*c+j);
            int at = (i*c+j)*CHANNELS;
            int top = ((itop)*c+j)*CHANNELS;
            int bottom = ((ibottom)*c+j)*CHANNELS;
            int left = (i*c+jleft)*CHANNELS;
            int right = (i*c+jright)*CHANNELS;
            _weights_inner(im,imC,at,real_at,top,left,right,bottom,Wn);
        }
    }
}

FORCE_INLINE void _weights_block(
        double *tmp_weights,
        int rSTART, int rEND, int cSTART, int cEND,
        double *im, double **imC, double *Wn,
        uint32_t r, uint32_t c){
    // copy
    for(int i = rSTART; i < rEND; i++){
        int small_first = ((i-rSTART)*(cEND-cSTART))*CHANNELS;
        int first = (i*c+cSTART)*CHANNELS;
        int len = (cEND-cSTART)*CHANNELS*sizeof(double);
        memcpy(&(tmp_weights[small_first]),&(im[first]),len);
    }

    //center rows
    for(int i = rSTART+1; i < rEND-1; i++){
        int ti = i-rSTART;
        int itop = ti-1;
        int ibottom = ti+1;
        //center cols
        for(int j = cSTART+1; j < cEND-1; j++){
            int tj = j-cSTART;
            int jleft = tj-1;
            int jright = tj+1;

            int real_at = (i*c+j);

            int tc = cEND-cSTART;
            int at = (ti*tc+tj)*CHANNELS;

            int top = ((itop)*tc+tj)*CHANNELS;
            int bottom = ((ibottom)*tc+tj)*CHANNELS;
            int left = (ti*tc+jleft)*CHANNELS;
            int right = (ti*tc+jright)*CHANNELS;

            _weights_inner(tmp_weights,imC,at,real_at,top,left,right,bottom,Wn);

        }
    }
}

FORCE_INLINE void _weights_block_unroll2x(
        double *tmp_weights,
        int rSTART, int rEND, int cSTART, int cEND,
        double *im, double **imC, double *Wn,
        uint32_t r, uint32_t c){
    // copy
    for(int i = rSTART; i < rEND; i++){
        int small_first = ((i-rSTART)*(cEND-cSTART))*CHANNELS;
        int first = (i*c+cSTART)*CHANNELS;
        int len = (cEND-cSTART)*CHANNELS*sizeof(double);
        memcpy(&(tmp_weights[small_first]),&(im[first]),len);
    }

    int tc = cEND-cSTART;
    double rexp1 = 0;
    double gexp1 = 0;
    double bexp1 = 0;
    double rexp2 = 0;
    double gexp2 = 0;
    double bexp2 = 0;
    __attribute__ ((aligned (16)))
    __attribute__ ((aligned (16))) double target[2];
    __attribute__ ((aligned (16))) double targetr[2];
    __attribute__ ((aligned (16))) double targetg[2];
    __attribute__ ((aligned (16))) double targetb[2];

    __m256d greyweights = _mm256_set_pd (0, 0.1140, 0.5870, 0.2989); //rgb0
    __m128d three = _mm_set1_pd(3.0);
    __m128d onehalf = _mm_set1_pd(0.5);
    __m128d eps = _mm_set1_pd(1.0E-12);
    __m128d minustwelvehalf = _mm_set1_pd(-12.5);
    __m128d signbit_mask128d = _mm_set1_pd(-0.); // -0. = 1 << 63
    __m256d signbit_mask256d = _mm256_set1_pd(-0.); // -0. = 1 << 63

//    printf("_weights_block_unroll2x: rSTART: %d rEND: %d\n",rSTART,rEND);

    //center rows
    for(int ti = 1; ti < (rEND-rSTART)-1; ti+=WEIGHTS_UNROLL_NB_R){
        int i = ti+rSTART;
        int itop = ti-1;
        int ibottom = ti+1;
        //center cols
        for(int tj = 1; tj+(WEIGHTS_UNROLL_NB_C-1) < (cEND-cSTART)-1; tj+=WEIGHTS_UNROLL_NB_C){
            int j = tj+cSTART;
            int real_at = (i*c+j);
            int at = (ti*tc+tj)*CHANNELS;

            int top = ((itop)*tc+tj)*CHANNELS;
            int bottom = ((ibottom)*tc+tj)*CHANNELS;

            __m256d top1 = _mm256_loadu_pd (&(tmp_weights[top]));
            __m256d top2 = _mm256_loadu_pd (&(tmp_weights[top+3]));

            __m256d bot1 = _mm256_loadu_pd (&(tmp_weights[bottom]));
            __m256d bot2 = _mm256_loadu_pd (&(tmp_weights[bottom+3]));

            __m256d left0 = _mm256_loadu_pd (&(tmp_weights[at-3]));
            __m256d cent1 = _mm256_loadu_pd (&(tmp_weights[at]));
            __m256d cent2 = _mm256_loadu_pd (&(tmp_weights[at+3]));
            __m256d right3 = _mm256_loadu_pd (&(tmp_weights[at+6]));


            //for streaming
            __m128d rg1 = _mm256_extractf128_pd(cent1,0); // lo
            __m128d rg2 = _mm256_extractf128_pd(cent2,0); // lo
            __m128d bX1 = _mm256_extractf128_pd(cent1,1); // hi
            __m128d bX2 = _mm256_extractf128_pd(cent2,1); // hi

            __m128d rr = _mm_unpacklo_pd(rg1,rg2);
            __m128d gg = _mm_unpackhi_pd(rg1,rg2);
            __m128d bb = _mm_unpacklo_pd(bX1,bX2);

            //grey conversion
            __m256d g_top1 = _mm256_mul_pd(top1,greyweights);
            __m256d g_top2 = _mm256_mul_pd(top2,greyweights);
            __m256d g_bot1 = _mm256_mul_pd(bot1,greyweights);
            __m256d g_bot2 = _mm256_mul_pd(bot2,greyweights);
            __m256d g_left0 = _mm256_mul_pd(left0,greyweights);
            __m256d g_cent1 = _mm256_mul_pd(cent1,greyweights);
            __m256d g_cent2 = _mm256_mul_pd(cent2,greyweights);
            __m256d g_right3 = _mm256_mul_pd(right3,greyweights);

            __m256d s_t1b1 = _mm256_hadd_pd(g_top1,g_bot1);
            __m256d s_t2b2 = _mm256_hadd_pd(g_top2,g_bot2);
            __m256d s_l0c2 = _mm256_hadd_pd(g_left0,g_cent2);
            __m256d s_c1r3 = _mm256_hadd_pd(g_cent1,g_right3);

            __m256d v_t1b1t2b2 = _mm256_hadd_pd(s_t1b1,s_t2b2);
            __m256d v_l0c2c1r3 = _mm256_hadd_pd(s_l0c2,s_c1r3);

            //stddev
            __m128d sum_rg = _mm_add_pd(rr,gg);
            __m128d sum_rgb = _mm_add_pd(sum_rg,bb);
            __m128d mu12 = _mm_div_pd(sum_rgb,three);

            __m128d rmu = _mm_sub_pd(rr,mu12);
            __m128d gmu = _mm_sub_pd(gg,mu12);
            __m128d bmu = _mm_sub_pd(bb,mu12);

            __m128d rz = _mm_sub_pd(rr,onehalf);
            __m128d gz = _mm_sub_pd(gg,onehalf);
            __m128d bz = _mm_sub_pd(bb,onehalf);

            __m128d rzrz = _mm_mul_pd(rz,rz);
            __m128d gzgz = _mm_mul_pd(gz,gz);
            __m128d bzbz = _mm_mul_pd(bz,bz);

            __m128d rzrzt = _mm_mul_pd(rzrz,minustwelvehalf);
            __m128d gzgzt = _mm_mul_pd(gzgz,minustwelvehalf);
            __m128d bzbzt = _mm_mul_pd(bzbz,minustwelvehalf);

            _mm_storel_pd(&rexp1, rzrzt);
            _mm_storel_pd(&gexp1, gzgzt);
            _mm_storel_pd(&bexp1, bzbzt);
            _mm_storeh_pd(&rexp2, rzrzt);
            _mm_storeh_pd(&gexp2, gzgzt);
            _mm_storeh_pd(&bexp2, bzbzt);

            __m128d re = _mm_set_pd(exp(rexp2), exp(rexp1));
            __m128d ge = _mm_set_pd(exp(gexp2), exp(gexp1));
            __m128d be = _mm_set_pd(exp(bexp2), exp(bexp1));

            __m128d rmurmu = _mm_mul_pd(rmu,rmu);
            __m128d gmugmu = _mm_mul_pd(gmu,gmu);
            __m128d bmubmu = _mm_mul_pd(bmu,bmu);

            __m128d sum_mu_rg = _mm_add_pd(rmurmu,gmugmu);
            __m128d sum_mu_rgb = _mm_add_pd(sum_mu_rg,bmubmu);
            __m128d dev12 = _mm_div_pd(sum_mu_rgb,three);
            __m128d stddev12 = _mm_sqrt_pd(dev12);

            __m128d abs_stddev12 = _mm_andnot_pd(signbit_mask128d, stddev12); // !sign_mask & x

            __m128d rege12 = _mm_mul_pd(re,ge);
            __m128d regebe12 = _mm_mul_pd(rege12,be);
            __m128d abs_regebe12 = _mm_andnot_pd(signbit_mask128d, regebe12); // !sign_mask & x

            __m128d tmp1 = _mm_mul_pd(abs_stddev12, abs_regebe12);

            //convolution step
//            __m256d v_t1b1t2b2
//            __m256d v_l0c2c1r3
//            r1 = (t1 + b1) + -4.0*c1 + (l0 + c2);
//            r2 = (t2 + b2) + -4.0*c2 + (c1 + r3);

            __m256d v_c2c2c1c1 = _mm256_shuffle_pd(v_l0c2c1r3,v_l0c2c1r3, _MM_SHUFFLE(2,2,1,1));
            __m256d neg_c2c2c1c1 = _mm256_xor_pd(signbit_mask256d, v_c2c2c1c1);
            __m256d v_minus4c12 = _mm256_hadd_pd(neg_c2c2c1c1, neg_c2c2c1c1);
            //  [-2*c2,-2*c1, -2*c2,-2*c1]
            __m256d v_2c2_2c2_2c1_2c1 = _mm256_unpacklo_pd(v_minus4c12,v_minus4c12);


            __m256d v_tb12_lcr12 = _mm256_hadd_pd(v_t1b1t2b2, v_l0c2c1r3);
            // [tb1, tb2, lcr1, lcr2]

            __m256d v_tb1_tb2 = _mm256_unpacklo_pd(v_tb12_lcr12, v_tb12_lcr12);
            // [tb1 tb1 tb2 tb2]
            __m256d v_lcr1_lcr2 = _mm256_unpackhi_pd(v_tb12_lcr12, v_tb12_lcr12);
            // [lcr1 lcr1 lcr2 lcr2]
            __m256d v_n1_n1_n2_n2 = _mm256_add_pd(v_tb1_tb2, v_lcr1_lcr2);
            //neighbors
            __m256d v_convolution = _mm256_hadd_pd(v_2c2_2c2_2c1_2c1, v_2c2_2c2_2c1_2c1 );
            // [c2 c1 c2 c1]
            __m256d v_rX_r1_r2_rX = _mm256_add_pd(v_n1_n1_n2_n2, v_convolution);
            // [ _ r1 r2 _ }

            __m128d Xr1 = _mm256_extractf128_pd (v_rX_r1_r2_rX, 0);
            __m128d r2X = _mm256_extractf128_pd (v_rX_r1_r2_rX, 1);

            __m128d r1r2 = _mm_shuffle_pd (Xr1, r2X, 1); //01 mask
            __m128d abs_r1r2 = _mm_andnot_pd(signbit_mask128d, r1r2); // !sign_mask & x
            __m128d tmp2 = _mm_mul_pd(tmp1,abs_r1r2);
            __m128d tmp3 = _mm_add_pd(tmp2,eps);


            //store weight
            //TODO streaming
            _mm_store_pd(&(target[0]),tmp3);

            //store individual channels
            //TODO streaming
            _mm_store_pd (&(targetr[0]), rr);
            _mm_store_pd (&(targetg[0]), gg);
            _mm_store_pd (&(targetb[0]), bb);

            memcpy(&(imC[0][real_at]),&(targetr[0]),2*sizeof(double));
            memcpy(&(imC[1][real_at]),&(targetg[0]),2*sizeof(double));
            memcpy(&(imC[2][real_at]),&(targetb[0]),2*sizeof(double));
            memcpy(&(Wn[real_at]),&(target[0]),2*sizeof(double));

            COST_INC_ADD(18);
            COST_INC_MUL(17);
            COST_INC_DIV(2);
            COST_INC_ABS(3);
            COST_INC_SQRT(1);
            COST_INC_EXP(3);
            COST_INC_ADD(18);
            COST_INC_MUL(17);
            COST_INC_DIV(2);
            COST_INC_ABS(3);
            COST_INC_SQRT(1);
            COST_INC_EXP(3);
        }
    }
}


