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
