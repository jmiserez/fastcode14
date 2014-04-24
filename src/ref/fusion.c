#include<fusion.h>
#include"segments.h"

void ones(double *dst, size_t len);
void contrast(double *im, int w, int h, double *C);

double* exposure_fusion(double** I, int w, int h, int N, double m[3],
    void* _segments) {
    int i;
    size_t npixels = w * h;

    segments_t* mem = (segments_t*) _segments;

    double contrast_parm = m[0];
    double sat_parm = m[1];
    double wexp_parm = m[2];

    double* temp_matrix = mem->temp_matrix;
    double** weight_matrices = mem->weight_matrices;

    for( int i = 0; i < N; i++ ) {
        if( contrast_parm > 0.0 ) {

        }
    }
}

void ones(double *dst, size_t len){
    for(int i = 0; i < len; i++){
        dst[i] = (double)1.0;
    }
}

void contrast( double *C, int w, int h, double *im, segments_t* mem ){
    //laplacian filter
    double f[] = {
        0.0, 1.0, 0.0,
        1.0, -4.0, 1.0,
        0.0, 1.0, 0.0
    };

    double *mono_matrix = mem->mono_matrix;
    size_t npixels = w*h; //1 value/pixel
    zeros(C, npixels);
    // for each image, calculate contrast measure on grayscale version of the
    // image
    rgb2gray(im, npixels, mono_matrix);
    conv3x3_monochrome_replicate(mono_matrix,r,c,f,C);
}

/**
 * @brief convolution of a monochrome image with a 3x3 filter and border mode "replication"
 */
void conv3x3_monochrome_replicate( double* dst, double* im, int w, int h,
                                   double* f ) {
    for(int i = 1; i < h-1; i++) {
        for(int j = 1; j < w-1; j++){
            int iw0 = (i-1)*w, iw1 = (i)*w, iw2 = (i+1)*w;
            int j0 = j-1, j1 = j, j2 = j+1;
            dst[i*w+j] =
                    im[iw0+j0]*f[0] + im[iw0+j1]*f[1] + im[iw0+j2]*f[1] +
                    im[iw1+j0]*f[2] + im[iw1+j1]*f[3] + im[iw1+j2]*f[4] +
                    im[iw2+j0]*f[5] + im[iw1+j1]*f[6] + im[iw2+j3]*f[7];
        }
    }
    //edges
    for(int i = 1; i < h-1; i++){
        int j = 0;
        int iw0 = (i-1)*w, iw1 = (i)*w, iw2 = (i+1)*w;

        dst[i*w+j] =
                im[iw0+(j)]*f[0] + im[iw0+(j)]*f[1] + im[iw0+(j+1)]*f[1] +
                im[iw1+(j)]*f[2] + im[iw1+(j)]*f[3] + im[iw1+(j+1)]*f[4] +
                im[iw2+(j)]*f[5] + im[iw2+(j)]*f[6] + im[iw2+(j+1)]*f[7];
        j = w-1;
        dst[i*w+j] =
                im[iw0+(j-1)]*f[0] + im[iw0+(j)]*f[1] + im[iw0+(j)]*f[1] +
                im[iw1+(j-1)]*f[2] + im[iw1+(j)]*f[3] + im[iw1+(j)]*f[4] +
                im[iw2+(j-1)]*f[5] + im[iw2+(j)]*f[6] + im[iw2+(j)]*f[7];
    }
    for(int j = 1; j < w-1; j++){
        int i = 0;
        dst[i*w+j] =
                im[0+(j-1)]*f[0] + im[0+(j)]*f[1] + im[0+(j+1)]*f[1] +
                im[0+(j-1)]*f[2] + im[0+(j)]*f[3] + im[0+(j+1)]*f[4] +
                im[w+(j-1)]*f[5] + im[w+(j)]*f[6] + im[w+(j+1)]*f[7];
        i = h-1;
        int w0 = (h-2)*w, w1 = (h-1)*w;
        dst[i*w+j] =
                im[w0+(j-1)]*f[0] + im[w0+(j)]*f[1] + im[w0+(j+1)]*f[1] +
                im[w1+(j-1)]*f[2] + im[w1+(j)]*f[3] + im[w1+(j+1)]*f[4] +
                im[w1+(j-1)]*f[5] + im[w1+(j)]*f[6] + im[w1+(j+1)]*f[7];
    }
    //corners
    dst[0] =
            im[0]*f[0] + im[0]*f[1] + im[0  ]*f[1] +
            im[0]*f[2] + im[0]*f[3] + im[0  ]*f[4] +
            im[w]*f[5] + im[w]*f[6] + im[w+1]*f[7];
    dst[w-1] =
            im[w-2  ]*f[0] + im[w-1  ]*f[1] + im[w-1  ]*f[1] +
            im[w-2  ]*f[2] + im[w-1  ]*f[3] + im[w-1  ]*f[4] +
            im[2*w-2]*f[5] + im[2*w-1]*f[6] + im[2*w-1]*f[7];
    dst[(h-1)*w] =
            im[(h-2)*w]*f[0] + im[(h-2)*w]*f[1] + im[(h-2)*w+1]*f[1] +
            im[(h-1)*w]*f[2] + im[(h-1)*w]*f[3] + im[(h-1)*w+1]*f[4] +
            im[(h-1)*w]*f[5] + im[(h-1)*w]*f[6] + im[(h-1)*w+1]*f[7];
    int w0 = (h-1)*w-2, w1 = (h-1)*w-1;
    dst[h*w-1] =
            im[w0   ]*f[0] + im[w1   ]*f[1] + im[w1   ]*f[1] +
            im[h*w-2]*f[2] + im[h*w-1]*f[3] + im[h*w-1]*f[4] +
            im[h*w-2]*f[5] + im[h*w-1]*f[6] + im[h*w-1]*f[7];

}
