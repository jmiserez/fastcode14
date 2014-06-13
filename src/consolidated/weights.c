#ifndef WEIGHTS_C
#define WEIGHTS_C

#include "weights.h"

FORCE_INLINE void _weights_pixel(double *im, int at3, int at, int top, int left, int right, int bottom, double *Wn, double sat_parm, double contrast_parm, double wexp_parm){
    double r = im[at3];
    double g = im[at3+1];
    double b = im[at3+2];
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
    double t2 = pow(fabs(t1),sat_parm);
    COST_INC_POW(1);
    COST_INC_ABS(1);
    double t3 = re*ge*be;
    COST_INC_MUL(2);
    double t4 = pow(fabs(t3),wexp_parm);
    COST_INC_POW(1);
    COST_INC_ABS(1);

    double t5 = t2 * t4;
    COST_INC_MUL(1);

    double t6 = -4.0*grey+grey1+grey2+grey3+grey4;
    COST_INC_MUL(1);
    COST_INC_ADD(4);

    double t7 = pow(fabs(t6),contrast_parm);
    COST_INC_POW(1);
    COST_INC_ABS(1);

    double t8 = t5 * t7;
    COST_INC_MUL(1);

    double t9 = t8 + 1.0E-12;
    COST_INC_ADD(1);

    Wn[at] = t9;
}

FORCE_INLINE void _weights_inner_row(double *im, int c, int i, int itop, int ibottom, double *Wn, double sat_parm, double contrast_parm, double wexp_parm){

    //center cols
    for(int j = 1; j < c-1; j++){
        double at = (i*c+j);
        double at3 = at*CHANNELS;
        int top = ((itop)*c+j)*CHANNELS;
        int bottom = ((ibottom)*c+j)*CHANNELS;
        int left = (i*c+j-1)*CHANNELS;
        int right = (i*c+j+1)*CHANNELS;
        _weights_pixel(im,at3,at,top,left,right,bottom,Wn,sat_parm,contrast_parm,wexp_parm);
    }
    //left cols
    int j = 0;
    double at = (i*c+j);
    double at3 = at*CHANNELS;
    int top = ((itop)*c+j)*CHANNELS;
    int bottom = ((ibottom)*c+j)*CHANNELS;
    int left = (i*c+j)*CHANNELS;
    int right = (i*c+j+1)*CHANNELS;
    _weights_pixel(im,at3,at,top,left,right,bottom,Wn,sat_parm,contrast_parm,wexp_parm);
    //right cols
    j = c-1;
    at = (i*c+j);
    at3 = at*CHANNELS;
    top = ((itop)*c+j)*CHANNELS;
    bottom = ((ibottom)*c+j)*CHANNELS;
    left = (i*c+j-1)*CHANNELS;
    right = (i*c+j)*CHANNELS;
    _weights_pixel(im,at3,at,top,left,right,bottom,Wn,sat_parm,contrast_parm,wexp_parm);
}

void weights(uint32_t nimages, uint32_t r, uint32_t c, double contrast_parm, double sat_parm, double wexp_parm,
             double **I, double **W){

    for (int n = 0; n < nimages; n++){
        double *im = I[n];
        double *Wn = W[n];

        //center rows
        for(int i = 1; i < r-1; i++){
            int itop = i-1;
            int ibottom = i+1;
            _weights_inner_row(im,c,i,itop,ibottom,Wn,sat_parm,contrast_parm,wexp_parm);
        }
        //top row
        int i = 0;
        int itop = 0;
        int ibottom = 1;
        _weights_inner_row(im,c,i,itop,ibottom,Wn,sat_parm,contrast_parm,wexp_parm);
        //bottom row
        i = r-1;
        itop = r-2;
        ibottom = r-1;
        _weights_inner_row(im,c,i,itop,ibottom,Wn,sat_parm,contrast_parm,wexp_parm);
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
