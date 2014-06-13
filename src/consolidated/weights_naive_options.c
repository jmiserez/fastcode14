#ifndef WEIGHTS_C
#define WEIGHTS_C

#include "weights_naive_options.h"

void weights(uint32_t nimages, uint32_t npixels, uint32_t r, uint32_t c, double contrast_parm, double sat_parm, double wexp_parm,
             double **I, double **W, double *tmp_weights, double *tmp2_weights){
    //for each image, calculate the weight maps
    for (int n = 0; n < nimages; n++){
        for(int i = 0; i < npixels; i++){
            W[n][i] = (double)1.0;
        }

        if(contrast_parm > 0){
            rgb2gray(I[n], npixels, tmp2_weights); //TODO optimize
            conv3x3_monochrome_replicate(tmp2_weights,r,c,tmp_weights);
            for(int i = 0; i < npixels; i++){
                tmp_weights[i] = pow(fabs(tmp_weights[i]),contrast_parm); COST_INC_POW(1); COST_INC_ABS(1);
            }
            for(int i = 0; i < npixels; i++){
                W[n][i] = W[n][i] * tmp_weights[i]; COST_INC_MUL(1);
            }
        }

        if(sat_parm > 0){
            //saturation is computed as the standard deviation of the color channels
            for(int i = 0; i < npixels; i++){
                double r = I[n][i*3];
                double g = I[n][i*3+1];
                double b = I[n][i*3+2];
                double mu = (r + g + b) / 3.0;
                COST_INC_ADD(2);
                COST_INC_DIV(1);
                double rmu = r-mu;
                double gmu = g-mu;
                double bmu = b-mu;
                COST_INC_ADD(3);
                tmp_weights[i] = sqrt((rmu*rmu + gmu*gmu + bmu*bmu)/3.0);
                COST_INC_SQRT(1);
                COST_INC_ADD(3);
                COST_INC_MUL(3);
                COST_INC_DIV(1);
            }
            for(int i = 0; i < npixels; i++){
                tmp_weights[i] = pow(fabs(tmp_weights[i]),sat_parm); COST_INC_POW(1); COST_INC_ABS(1);
            }
            for(int i = 0; i < npixels; i++){
                W[n][i] = W[n][i] * tmp_weights[i]; COST_INC_MUL(1);
            }
        }

        if(wexp_parm > 0){
            for(int i = 0; i < npixels; i++){
                double r = I[n][i*3];
                double g = I[n][i*3+1];
                double b = I[n][i*3+2];
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
                tmp_weights[i] = r*g*b;
                COST_INC_MUL(2);
            }
            for(int i = 0; i < npixels; i++){
                tmp_weights[i] = pow(fabs(tmp_weights[i]),wexp_parm); COST_INC_POW(1); COST_INC_ABS(1);
            }
            for(int i = 0; i < npixels; i++){
                W[n][i] = W[n][i] * tmp_weights[i]; COST_INC_MUL(1);
            }
        }

        for(int i = 0; i < npixels; i++){
            W[n][i] = W[n][i] + 1.0E-12; COST_INC_ADD(1);
        }
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

/**
 * @brief Implementation of the MATLAB rgb2gray function
 *
 * See: http://www.mathworks.com/help/images/ref/rgb2gray.html
 *
 * @param rgb Input image
 * @param npixels Size of image in pixels
 * @param gray (out) Output image
 */
void rgb2gray(double *im, size_t npixels, double* dst){
    for(int i = 0; i < npixels; i++){
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        dst[i] = 0.2989 * r + 0.5870 * g + 0.1140 * b; //retain luminance, discard hue and saturation
        COST_INC_ADD(2);
        COST_INC_MUL(3);
    }
}

/**
 * @brief convolution of a monochrome image with a 3x3 filter and border mode "replication"
 */
void conv3x3_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* dst){
    for(int i = 1; i < r-1; i++){
        for(int j = 1; j < c-1; j++){
            dst[i*c+j] =
                    im[(i-1)*c+(j)] +
                    im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                    im[(i+1)*c+(j)];
             COST_INC_ADD(4);
             COST_INC_MUL(1);
        }
    }
    //edges
    for(int i = 1; i < r-1; i++){
        int j = 0;
        dst[i*c+j] =
                im[(i-1)*c+(j)] +
                im[(i)  *c+(j)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                im[(i+1)*c+(j)];
        COST_INC_ADD(4);
        COST_INC_MUL(1);
        j = c-1;
        dst[i*c+j] =
                im[(i-1)*c+(j)] +
                im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j)] +
                im[(i+1)*c+(j)];
        COST_INC_ADD(4);
        COST_INC_MUL(1);
    }
    for(int j = 1; j < c-1; j++){
        int i = 0;
        dst[i*c+j] =
                im[(i)  *c+(j)] +
                im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                im[(i+1)*c+(j)];
        COST_INC_ADD(4);
        COST_INC_MUL(1);
        i = r-1;
        dst[i*c+j] =
                im[(i-1)*c+(j)] +
                im[(i)  *c+(j-1)] + im[(i)  *c+(j)]*-4.0 + im[(i)  *c+(j+1)] +
                im[(i)  *c+(j)];
        COST_INC_ADD(4);
        COST_INC_MUL(1);
    }
    //corners
    //top left
    int i = 0;
    int j = 0;
    dst[i*c+j] =
            im[(i  )*c+(j)] +
            im[(i  )*c+(j)] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j+1)] +
            im[(i+1)*c+(j)];
    COST_INC_ADD(4);
    COST_INC_MUL(1);
    //top right
    i = 0;
    j = c-1;
    dst[i*c+j] =
            im[(i  )*c+(j)] +
            im[(i  )*c+(j-1)] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j  )] +
            im[(i+1)*c+(j)];
    COST_INC_ADD(4);
    COST_INC_MUL(1);
    //bottom left
    i = r-1;
    j = 0;
    dst[i*c+j] =
            im[(i-1)*c+(j)] +
            im[(i  )*c+(j  )] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j+1)] +
            im[(i  )*c+(j)];
    COST_INC_ADD(4);
    COST_INC_MUL(1);
    //bottom right
    i = r-1;
    j = c-1;
    dst[i*c+j] =
            im[(i-1)*c+(j)] +
            im[(i  )*c+(j-1)] + im[(i  )*c+(j)]*-4.0 + im[(i  )*c+(j  )] +
            im[(i  )*c+(j)];
    COST_INC_ADD(4);
    COST_INC_MUL(1);

}

#endif
