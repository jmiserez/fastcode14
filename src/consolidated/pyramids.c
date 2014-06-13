#include "pyramids.h"

void pyramids(uint32_t nimages, uint32_t nlev, uint32_t r, uint32_t c, double **I, double **W,
              double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize,
              double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
              double ***pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c){
    //multiresolution blending
    for (int n = 0; n < nimages; n++){
        //construct 1-channel gaussian pyramid from weights
        gaussian_pyramid(W[n],r,c,1,nlev,tmp_halfsize,pyrW[n],pyrW_r[n],pyrW_c[n]);
    }
    for (int n = 0; n < nimages; n++){
        //construct 3-channel laplacian pyramid from images
        laplacian_pyramid(I[n],r,c,nlev,tmp_halfsize,tmp_quartsize,tmp2_quartsize,pyrI[n],pyrI_r[n],pyrI_c[n]);
    }
}

void blend(uint32_t nimages, uint32_t nlev,
           double **pyr, uint32_t *pyr_r, uint32_t *pyr_c,
           double ***pyrW, uint32_t **pyrW_r, uint32_t **pyrW_c,
           double ***pyrI, uint32_t **pyrI_r, uint32_t **pyrI_c){
    //weighted blend
    if(0 < nimages){
        int n = 0;
        for(int v = 0; v < nlev; v++){
            for(int i = 0; i < pyrI_r[n][v]; i++){
                for(int j = 0; j < pyrI_c[n][v]; j++){
                    int p_idx3 = 3*(i*pyr_c[v]+j);
                    int pI_idx3 = 3*(i*pyrI_c[n][v]+j);
                    int pW_idx = i*pyrI_c[n][v]+j;

                    pyr[v][p_idx3] = pyrW[n][v][pW_idx] * pyrI[n][v][pI_idx3];
                    pyr[v][p_idx3+1] = pyrW[n][v][pW_idx] * pyrI[n][v][pI_idx3+1];
                    pyr[v][p_idx3+2] = pyrW[n][v][pW_idx] * pyrI[n][v][pI_idx3+2];
                    COST_INC_MUL(3);
                }
            }
        }
    }
    for (int n = 1; n < nimages; n++){
        for(int v = 0; v < nlev; v++){
            for(int i = 0; i < pyrI_r[n][v]; i++){
                for(int j = 0; j < pyrI_c[n][v]; j++){
                    int p_idx3 = 3*(i*pyr_c[v]+j);
                    int pI_idx3 = 3*(i*pyrI_c[n][v]+j);
                    int pW_idx = i*pyrI_c[n][v]+j;

                    pyr[v][p_idx3+] += pyrW[n][v][pW_idx] * pyrI[n][v][pI_idx3];
                    pyr[v][p_idx3+1] += pyrW[n][v][pW_idx] * pyrI[n][v][pI_idx3+1];
                    pyr[v][p_idx3+2] += pyrW[n][v][pW_idx] * pyrI[n][v][pI_idx3+2];
                    COST_INC_ADD(3);
                    COST_INC_MUL(3);
                }
            }
        }
    }
}

void reconstruct_laplacian_pyramid(uint32_t nlev, double *tmp_fullsize, double *tmp2_fullsize, double *tmp_halfsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c, uint32_t r, uint32_t c, double *dst){
    if (nlev-2 >= 0){
        int v = nlev-2;
        upsample(pyr[v+1],pyr_r[v+1],pyr_c[v+1],3,pyr_r[v],pyr_c[v],tmp_halfsize,tmp2_fullsize);
        for(int i = 0; i < pyr_r[v]*pyr_c[v]*3; i++){
            dst[i] = pyr[v][i] + tmp2_fullsize[i]; COST_INC_ADD(1);
        }
    }
    for (int v = nlev-3; v >= 0; v--){
        upsample(dst,pyr_r[v+1],pyr_c[v+1],3,pyr_r[v],pyr_c[v],tmp_halfsize,tmp2_fullsize);
        for(int i = 0; i < pyr_r[v]*pyr_c[v]*3; i++){
            dst[i] = pyr[v][i] + tmp2_fullsize[i]; COST_INC_ADD(1);
        }
    }
}

void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    pyr[0] = im;
    if(1 < nlev){
        int v = 1;
        downsample(pyr[0],pyr_r[v-1],pyr_c[v-1],3,tmp_halfsize,pyr_r[v],pyr_c[v],pyr[v]);
    }
    for(int v = 2; v < nlev; v++){
        //downsample image and store into level
        downsample(pyr[v-1],pyr_r[v-1],pyr_c[v-1],3,tmp_halfsize,pyr_r[v],pyr_c[v],pyr[v]);
    }
}

void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double *tmp_halfsize, double *tmp_quartsize, double *tmp2_quartsize, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    uint32_t S_r = r;
    uint32_t S_c = c;
    uint32_t T_r = r;
    uint32_t T_c = c;

    double *tmp = NULL;

    if(0 < nlev-1){
        int v = 0;
        S_r = pyr_r[v+1];
        S_c = pyr_c[v+1];
        downsample(im,T_r,T_c,3,tmp_halfsize,S_r,S_c,tmp2_quartsize);
        upsample(tmp2_quartsize,S_r,S_c,3,pyr_r[v],pyr_c[v],tmp_halfsize,pyr[v]);
        for(int i = 0; i < T_r*T_c*3; i++){
            pyr[v][i] = im[i] - pyr[v][i]; COST_INC_ADD(1);
        }
        T_r = S_r;
        T_c = S_c;
        double *tmp = tmp_quartsize;
        tmp_quartsize = tmp2_quartsize;
        tmp2_quartsize = tmp;
    }
    for(int v = 1; v < nlev-1; v++){
        S_r = pyr_r[v+1];
        S_c = pyr_c[v+1];
        downsample(tmp_quartsize,T_r,T_c,3,tmp_halfsize,S_r,S_c,tmp2_quartsize);
        upsample(tmp2_quartsize,S_r,S_c,3,pyr_r[v],pyr_c[v],tmp_halfsize,pyr[v]);
        for(int i = 0; i < T_r*T_c*3; i++){
            pyr[v][i] = tmp_quartsize[i] - pyr[v][i]; COST_INC_ADD(1);
        }
        T_r = S_r;
        T_c = S_c;
        tmp = tmp_quartsize;
        tmp_quartsize = tmp2_quartsize;
        tmp2_quartsize = tmp;
    }
    //memcpy(dst,src,src_len*sizeof(double));
    for(int i = 0; i < T_r*T_c*3; i++){
        pyr[nlev-1][i] = tmp_quartsize[i];
    }
}

/**
 * @brief convolution of a multi-channel image with a separable 5x5 filter and border mode "symmetric" combined with downsampling (1/4th of the pixels)
 */
void downsample(double *im, uint32_t r, uint32_t c, double *tmp_halfsize, uint32_t down_r, uint32_t down_c, double *dst){
    //r is height (vertical), c is width (horizontal)
    //horizontal filter
    int c2 = c/2+((c-2) % 2); //tmp_halfsize is only half the size
    for(int i = 0; i < r; i++){
        for(int j = 2; j < c-2; j+=2){ //every 2nd column
                int tidx3 = 3*(i*c2+(j/2));
                int i1_idx3 = 3*((i  )*c+(j-2));
                int i2_idx3 = 3*((i  )*c+(j-1));
                int i3_idx3 = 3*((i  )*c+(j ));
                int i4_idx3 = 3*((i  )*c+(j+1));
                int i5_idx3 = 3*((i  )*c+(j+2));

                tmp_halfsize[tidx3] =
                        im[i1_idx3]*.0625 +
                        im[i2_idx3]*.25 +
                        im[i3_idx3]*.375 +
                        im[i4_idx3]*.25 +
                        im[i5_idx3]*.0625;
                tmp_halfsize[tidx3+1] =
                        im[i1_idx3+1]*.0625 +
                        im[i2_idx3+1]*.25 +
                        im[i3_idx3+1]*.375 +
                        im[i4_idx3+1]*.25 +
                        im[i5_idx3+1]*.0625;
                tmp_halfsize[tidx3+2] =
                        im[i1_idx3+2]*.0625 +
                        im[i2_idx3+2]*.25 +
                        im[i3_idx3+2]*.375 +
                        im[i4_idx3+2]*.25 +
                        im[i5_idx3+2]*.0625;
                COST_INC_ADD(12);
                COST_INC_MUL(15);
        }
        //left edge
        int j = 0; // 1 0 [0 1 2 ... ]
        for(int k = 0; k < 3; k++){
            tmp_halfsize[(i*c2+(j/2))*3+k] =
                    im[((i  )*c+(j+1))*3+k]*.3125 +
                    im[((i  )*c+(j  ))*3+k]*.625 +
                    im[((i  )*c+(j+2))*3+k]*.0625;
            COST_INC_ADD(2);
            COST_INC_MUL(3);
        }
        //right edge
        if((c-2) % 2 == 0){
            j = c-2; // [ ... -2 -1 0 1] 1
            for(int k = 0; k < 3; k++){
                tmp_halfsize[(i*c2+(j/2))*3+k] =
                        im[((i  )*c+(j-2))*3+k]*.0625 +
                        im[((i  )*c+(j-1))*3+k]*.25 +
                        im[((i  )*c+(j  ))*3+k]*.375 +
                        im[((i  )*c+(j+1))*3+k]*.3125;
                COST_INC_ADD(3);
                COST_INC_MUL(4);
            }
        }else{
            j = c-1; // [ ... -2 -1 0] 0 -1
            for(int k = 0; k < 3; k++){
                tmp_halfsize[(i*c2+(j/2))*3+k] =
                        im[((i  )*c+(j-2))*3+k]*.0625 +
                        im[((i  )*c+(j-1))*3+k]*.3125 +
                        im[((i  )*c+(j  ))*3+k]*.625;
                COST_INC_ADD(2);
                COST_INC_MUL(3);
            }
        }
    }
    //vertical filter
    c = c2;
    for(int j = 0; j < c; j++){ //every column in tmp_halfsize = every 2nd column in im
        for(int i = 2; i < r-2; i+=2){ //every 2nd row in tmp_halfsize
            int i2 = i/2+((i-2) % 2);

            int tidx3 = 3*(i2*c+j);
            int i1_idx3 = 3*((i-2)*c+(j));
            int i2_idx3 = 3*((i-2)*c+(j));
            int i3_idx3 = 3*((i  )*c+(j));
            int i4_idx3 = 3*((i+1)*c+(j));
            int i5_idx3 = 3*((i+2)*c+(j));

            dst[tidx3+] =
                    tmp_halfsize[i1_idx3+]*.0625 +
                    tmp_halfsize[i2_idx3+]*.25 +
                    tmp_halfsize[i3_idx3+]*.375 +
                    tmp_halfsize[i4_idx3+]*.25 +
                    tmp_halfsize[i5_idx3+]*.0625;
            dst[tidx3+2] =
                    tmp_halfsize[i1_idx3+1]*.0625 +
                    tmp_halfsize[i2_idx3+1]*.25 +
                    tmp_halfsize[i3_idx3+1]*.375 +
                    tmp_halfsize[i4_idx3+1]*.25 +
                    tmp_halfsize[i5_idx3+1]*.0625;
            dst[tidx3+2] =
                    tmp_halfsize[i1_idx3+2]*.0625 +
                    tmp_halfsize[i2_idx3+2]*.25 +
                    tmp_halfsize[i3_idx3+2]*.375 +
                    tmp_halfsize[i4_idx3+2]*.25 +
                    tmp_halfsize[i5_idx3+2]*.0625;
            COST_INC_ADD(12);
            COST_INC_MUL(15);
        }
        //top edge
        int i = 0; // 1 0 [0 1 2 ... ]
        int i2 = i/2+((i-2) % 2);
        for(int k = 0; k < 3; k++){
            dst[(i2*c+j)*3+k] =
                    tmp_halfsize[((i+1)*c+(j  ))*3+k]*.3125 +
                    tmp_halfsize[((i  )*c+(j  ))*3+k]*.625 +
                    tmp_halfsize[((i+2)*c+(j  ))*3+k]*.0625;
            COST_INC_ADD(2);
            COST_INC_MUL(3);
        }
        //bottom edge
        if((r-2) % 2 == 0){
            i = r-2; // [ ... -2 -1 0 1] 1
            i2 = i/2+((i-2) % 2);
            for(int k = 0; k < 3; k++){
                dst[(i2*c+j)*3+k] =
                        tmp_halfsize[((i-2)*c+(j  ))*3+k]*.0625 +
                        tmp_halfsize[((i-1)*c+(j  ))*3+k]*.25 +
                        tmp_halfsize[((i  )*c+(j  ))*3+k]*.375 +
                        tmp_halfsize[((i+1)*c+(j  ))*3+k]*.3125;
                COST_INC_ADD(3);
                COST_INC_MUL(4);
            }
        }else{
            i = r-1; // [ ... -2 -1 0] 0 -1
            i2 = i/2+((i-2) % 2);
            for(int k = 0; k < 3; k++){
                dst[(i2*c+j)*3+k] =
                        tmp_halfsize[((i-2)*c+(j  ))*3+k]*.0625 +
                        tmp_halfsize[((i-1)*c+(j  ))*3+k]*.3125 +
                        tmp_halfsize[((i  )*c+(j  ))*3+k]*.625;
                COST_INC_ADD(2);
                COST_INC_MUL(3);
            }
        }
    }
}

void upsample(double *im, uint32_t r, uint32_t c, uint32_t 3, uint32_t up_r, uint32_t up_c, double *tmp_halfsize, double *dst){
    // x
    for(int i = 0; i < r; i++){ //every 2nd line
        for(int j = 2; j < up_c-2; j++){
            for(int k = 0; k < 3; k++){
                tmp_halfsize[(i*up_c+j)*3+k] =
                        0.25*im[((i  )*c+(j/2-1))*3+k] +
                        1.5*im[((i )*c+(j/2  ))*3+k] +
                        0.25*im[((i  )*c+(j/2+1))*3+k];
                COST_INC_ADD(2);
                COST_INC_MUL(3);
            }
            j++;
            for(int k = 0; k < 3; k++){
                tmp_halfsize[(i*up_c+j)*3+k] =
                        im[((i  )*c+(j/2))*3+k] +
                        im[((i  )*c+(j/2+1))*3+k];
                COST_INC_ADD(1);
            }
        }
        //left edge
        int j = 0; // 0 0 [0 1 2 ... ]
        for(int k = 0; k < 3; k++){
            tmp_halfsize[(i*up_c+j)*3+k] =
                    1.75*im[((i  )*c+(j/2  ))*3+k] +
                    0.25*im[((i  )*c+(j/2+1))*3+k];
            COST_INC_ADD(1);
            COST_INC_MUL(2);
        }
        j = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < 3; k++){
            tmp_halfsize[(i*up_c+j)*3+k] =
                    im[((i  )*c+(j/2))*3+k] +
                    im[((i  )*c+(j/2+1))*3+k];
            COST_INC_ADD(1);
        }
        if(up_c % 2 == 0){
            //right edge
            j = up_c-2; // [ ... -2 -1 0 1] 0
            for(int k = 0; k < 3; k++){
                tmp_halfsize[(i*up_c+j)*3+k] =
                        0.25*im[((i  )*c+(j/2-1))*3+k] +
                        1.75*im[((i  )*c+(j/2  ))*3+k];
                COST_INC_ADD(1);
                COST_INC_MUL(2);
            }
            j = up_c-1; // [ ... -2 -1 0] 0 0
            for(int k = 0; k < 3; k++){
                tmp_halfsize[(i*up_c+j)*3+k] =
                        2.0*im[((i  )*c+(j/2))*3+k];
                COST_INC_MUL(1);
            }
        } else {
            //right edge (remaining)
            j = up_c-1; // [ ... -2 -1 0] 0 0
            for(int k = 0; k < 3; k++){
                tmp_halfsize[(i*up_c+j)*3+k] =
                        0.25*im[((i  )*c+(j/2-1))*3+k] +
                        1.75*im[((i  )*c+(j/2  ))*3+k];
                COST_INC_ADD(1);
                COST_INC_MUL(2);
            }
        }
    }

    // y
    for(int j = 0; j < up_c; j++){ //all columns
        for(int i = 1; i < r-1; i++){
            for(int k = 0; k < 3; k++){
                dst[(2*i*up_c+j)*3+k] =
                        tmp_halfsize[((i-1)*up_c+(j  ))*3+k]*.0625 +
                        tmp_halfsize[((i  )*up_c+(j  ))*3+k]*.375 +
                        tmp_halfsize[((i+1)*up_c+(j  ))*3+k]*.0625;
                COST_INC_ADD(2);
                COST_INC_MUL(3);
            }
            for(int k = 0; k < 3; k++){
                dst[((2*i+1)*up_c+j)*3+k] =
                        tmp_halfsize[((i)*up_c+(j  ))*3+k]*.25 +
                        tmp_halfsize[((i+1)*up_c+(j  ))*3+k]*.25;
                COST_INC_ADD(1);
                COST_INC_MUL(2);
            }
        }
        //top edge
        int i = 0; // 0 0 [0 1 2 ... ]
        for(int k = 0; k < 3; k++){
            dst[(2*0*up_c+j)*3+k] =
                    tmp_halfsize[((i  )*up_c+(j  ))*3+k]*.4375 +
                    tmp_halfsize[((i+1)*up_c+(j  ))*3+k]*.0625;
            COST_INC_ADD(1);
            COST_INC_MUL(2);
        }
        for(int k = 0; k < 3; k++){
            dst[((2*0+1)*up_c+j)*3+k] =
                    tmp_halfsize[((i)*up_c+(j  ))*3+k]*.25 +
                    tmp_halfsize[((i+1)*up_c+(j  ))*3+k]*.25;
            COST_INC_ADD(1);
            COST_INC_MUL(2);
        }
        if(up_r % 2 == 0){
            //bottom edge
            i = r-1; // [ ... -2 -1 0 1] 0
            for(int k = 0; k < 3; k++){
                dst[((2*i)*up_c+j)*3+k] =
                        tmp_halfsize[((i-1)*up_c+(j  ))*3+k]*.0625 +
                        tmp_halfsize[((i  )*up_c+(j  ))*3+k]*.4375;
                COST_INC_ADD(1);
                COST_INC_MUL(2);
            }
            // [ ... -2 -1 0] 0 0
            for(int k = 0; k < 3; k++){
                dst[((2*i+1)*up_c+j)*3+k] =
                        tmp_halfsize[((i)*up_c+(j  ))*3+k]*.5;
                COST_INC_MUL(1);
            }
        }else{
            //bottom edge
            i = r-1; // [ ... -2 -1 0] 0 0
            for(int k = 0; k < 3; k++){
                dst[((2*i)*up_c+j)*3+k] =
                        tmp_halfsize[((i-1)*up_c+(j  ))*3+k]*.0625 +
                        tmp_halfsize[((i  )*up_c+(j  ))*3+k]*.4375;
                COST_INC_ADD(1);
                COST_INC_MUL(2);
            }
        }
    }
}
