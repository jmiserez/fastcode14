#ifndef WEIGHTS_C
#define WEIGHTS_C

#include "weights_store_grey.h"

void weights(uint32_t nimages, uint32_t npixels, uint32_t r, uint32_t c,
             double **I, double **W, double *tmp_weights, double *tmp2_weights){
    //for each image, calculate the weight maps
    for (int n = 0; n < nimages; n++){
        for(int i = 0; i < npixels; i++){
            //saturation is computed as the standard deviation of the color channels
            double r = I[n][i*3];
            double g = I[n][i*3+1];
            double b = I[n][i*3+2];
            COST_INC_LOAD(3);
            double grey = 0.2989 * r + 0.5870 * g + 0.1140 * b; //retain luminance, discard hue and saturation
            COST_INC_ADD(2);
            COST_INC_MUL(3);
            double mu = (r + g + b) / 3.0;
            COST_INC_ADD(2);
            COST_INC_DIV(1);
            double rmu = r-mu;
            double gmu = g-mu;
            double bmu = b-mu;
            COST_INC_ADD(3);
            double factor_sat = fabs(sqrt((rmu*rmu + gmu*gmu + bmu*bmu)/3.0));
            COST_INC_MUL(1);
            COST_INC_SQRT(1);
            COST_INC_ADD(3);
            COST_INC_MUL(3);
            COST_INC_DIV(1);
            COST_INC_ABS(1);
            //well exposedness
            double rz = r-0.5;
            double gz = g-0.5;
            double bz = b-0.5;
            COST_INC_ADD(3);
            double rzrz = rz*rz;
            double gzgz = gz*gz;
            double bzbz = bz*bz;
            COST_INC_MUL(3);
            r = exp(-12.5*rzrz);
            g = exp(-12.5*gzgz);
            b = exp(-12.5*bzbz);
            COST_INC_MUL(3);
            COST_INC_EXP(3);
            double factor_wexp = fabs(r*g*b);
            COST_INC_MUL(1);
            COST_INC_ABS(1);
            COST_INC_MUL(2);
            //store weights for sat and wexp in separate array
            tmp2_weights[i] = factor_sat * factor_wexp;
            COST_INC_STORE(1);
            //store greyscale image
            tmp_weights[i] = grey;
            COST_INC_STORE(1);
        }
        conv3x3_mult_monochrome_replicate(tmp_weights,r,c,tmp2_weights,W[n]);
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

/**
 * @brief convolution of a monochrome image with a 3x3 filter and border mode "replication"
 *        multiplies the result with the a value from another location
 */
void conv3x3_mult_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* factors, double* dst){
    for(int i = 1; i < r-1; i++){
        for(int j = 1; j < c-1; j++){
            dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
                    im[(i-1)*c+(j)] +
                    im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                    im[(i+1)*c+(j)]);
             COST_INC_ADD(4);
             COST_INC_MUL(1);
             COST_INC_ABS(1);
             COST_INC_LOAD(6);
             COST_INC_STORE(1);
        }
    }
    //edges
    for(int i = 1; i < r-1; i++){
        int j = 0;
        dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
                im[(i-1)*c+(j)] +
                im[(i)  *c+(j)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                im[(i+1)*c+(j)]);
        COST_INC_ADD(4);
        COST_INC_MUL(1);
        COST_INC_ABS(1);
        COST_INC_LOAD(6);
        COST_INC_STORE(1);
        j = c-1;
        dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
                im[(i-1)*c+(j)] +
                im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j)] +
                im[(i+1)*c+(j)]);
        COST_INC_ADD(4);
        COST_INC_MUL(1);
        COST_INC_ABS(1);
        COST_INC_LOAD(6);
        COST_INC_STORE(1);
    }
    for(int j = 1; j < c-1; j++){
        int i = 0;
        dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
                im[(i)  *c+(j)] +
                im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                im[(i+1)*c+(j)]);
        COST_INC_ADD(4);
        COST_INC_MUL(1);
        COST_INC_ABS(1);
        COST_INC_LOAD(6);
        COST_INC_STORE(1);
        i = r-1;
        dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
                im[(i-1)*c+(j)] +
                im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                im[(i)  *c+(j)]);
        COST_INC_ADD(4);
        COST_INC_MUL(1);
        COST_INC_ABS(1);
        COST_INC_LOAD(6);
        COST_INC_STORE(1);
    }
    //corners
    //top left
    int i = 0;
    int j = 0;
    dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
            im[(i  )*c+(j)] +
            im[(i  )*c+(j)] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j+1)] +
            im[(i+1)*c+(j)]);
    COST_INC_ADD(4);
    COST_INC_MUL(1);
    COST_INC_ABS(1);
    COST_INC_LOAD(6);
    COST_INC_STORE(1);
    //top right
    i = 0;
    j = c-1;
    dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
            im[(i  )*c+(j)] +
            im[(i  )*c+(j-1)] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j  )] +
            im[(i+1)*c+(j)]);
    COST_INC_ADD(4);
    COST_INC_MUL(1);
    COST_INC_ABS(1);
    COST_INC_LOAD(6);
    COST_INC_STORE(1);
    //bottom left
    i = r-1;
    j = 0;
    dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
            im[(i-1)*c+(j)] +
            im[(i  )*c+(j  )] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j+1)] +
            im[(i  )*c+(j)]);
    COST_INC_ADD(4);
    COST_INC_MUL(1);
    COST_INC_ABS(1);
    COST_INC_LOAD(6);
    COST_INC_STORE(1);
    //bottom right
    i = r-1;
    j = c-1;
    dst[i*c+j] = 1.0E-12 + factors[i*c+j] * fabs(
            im[(i-1)*c+(j)] +
            im[(i  )*c+(j-1)] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j  )] +
            im[(i  )*c+(j)]);
    COST_INC_ADD(4);
    COST_INC_MUL(1);
    COST_INC_ABS(1);
    COST_INC_LOAD(6);
    COST_INC_STORE(1);
}

#endif
