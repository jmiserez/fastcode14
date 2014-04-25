#include<fusion.h>
#include"segments.h"


// weighting functions
void contrast( double *C, int w, int h, double *im, segments_t* mem );
void saturation( double *C, int w, int h, double *im, segments_t* mem );
void well_exposedness( double *C, int w, int h, double *im, segments_t* mem );

// MATLAB functions
void rgb2gray( double* dst, double *im, size_t npixels );

// helper functions
void matrix_set( double *dst, size_t len, double val );
void weight_matrix( double* weight_matrix, double* weights, int len,
                    double exp ) {
    int i;
    for( i = 0; i < len; i++ )
        weight_matrix[i] *= pow(weights[i], exp);
}

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
            contrast( temp_matrix, w, h, I[i], mem );
            weight_matrix( weight_matrices[i], temp_matrix, npixels,
                           contrast_parm );
        }

        if( sat_parm > 0.0 ) {
            saturation( temp_matrix, w, h, I[i], mem );
            weight_matrix( weight_matrices[i], temp_matrix, npixels,
                           sat_parm );
        }

        if( wexp_parm > 0.0 ) {
            well_exposedness( temp_matrix, w, h, I[i], mem );
            weight_matrix( weight_matrices[i], temp_matrix, npixels,
                           wexp_parm );
        }
    }
}

void saturation( double *C, int w, int h, double *im, segments_t* mem ) {
    int npixels = w * h;
    matrix_set( C, npixels, 0.0 );

    for( int i = 0; i < npixels; i++ ) {
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        double mu = (r + g + b) / 3.0;
        C[i] = sqrt(pow(r-mu,2) + pow(g-mu,2) + pow(b-mu,2)/3.0);
    }
}

void well_exposedness( double *C, int w, int h, double *im, segments_t* mem ) {
    int npixels = w * h, i;
    matrix_set( C, 0.0 );

    double sig = 0.2;
    for( i = 0; i < npixels; i++ ) {
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        r = exp(-0.5*pow(r - 0.5,2) / pow(sig,2));
        g = exp(-0.5*pow(g - 0.5,2) / pow(sig,2));
        b = exp(-0.5*pow(b - 0.5,2) / pow(sig,2));
        C[i] = r*g*b;
    }
}

void matrix_set( double *dst, size_t len, double val ){
    for(int i = 0; i < len; i++){
        dst[i] = val;
    }
}

void contrast( double *C, int w, int h, double *im, segments_t* mem ){
    size_t npixels = w*h; //1 value/pixel

    //laplacian filter
    double f[] = {
        0.0, 1.0, 0.0,
        1.0, -4.0, 1.0,
        0.0, 1.0, 0.0
    };

    double *mono_matrix = mem->mono_matrix;
    matrix_set( C, npixels, 0.0 );

    // for each image, calculate contrast measure on grayscale version of the
    // image
    rgb2gray( mono_matrix, im, npixels );
    conv3x3_monochrome_replicate( C, mono_matrix, w, h, f );
}



inline double apply_filter( double *im, int* idx, double* f, int len ) {
    int i;
    double res = 0.0;
    for( i = 0; i < len; i++ ) {
        res += im[idx[i]]*f[i];
    }
    return res;
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
            dst[i*w+j] = apply_filter( im, (int[]) {
                                      (iw0+j0), (iw0+j1), (iw0+j2),
                                      (iw1+j0), (iw1+j1), (iw1+j2),
                                      (iw2+j0), (iw1+j1), (iw2+j3)
                                  }, f, 9 );
        }
    }
    //edges
    for( int i = 1; i < h-1; i++ ) {
        int j = 0;
        int iw0 = (i-1)*w, iw1 = (i)*w, iw2 = (i+1)*w;

        dst[i*w+j] = apply_filter( im, (int[]) {
                                       iw0+j, iw0+j, iw0+j+1,
                                       iw1+j, iw1+j, iw1+j+1,
                                       iw2+j, iw2+j, iw2+j+1
                                   }, f, 9);
        j = w-1;
        dst[i*w+j] = apply_filter( im, (int[]) {
                                       iw0+(j-1), iw0+j, iw0+j,
                                       iw1+(j-1), iw1+j, iw1+j,
                                       iw2+(j-1), iw2+j, iw2+j
                                   }, f, 9 );
    }
    for(int j = 1; j < w-1; j++){
        int i = 0;
        dst[i*w+j] = apply_filter( im, (int[]) {
                                       j-1,   j,   j+1,
                                       j-1,   j,   j+1,
                                       w+j-1, w+j, w+j+1
                                   }, f, 9 );
        i = h-1;
        int w0 = (h-2)*w, w1 = (h-1)*w;
        dst[i*w+j] = apply_filter( im, (int[]) {
                                       w0+(j-1), w0+j, w0+j+1,
                                       w1+(j-1), w1+j, w1+j+1,
                                       w1+(j-1), w1+j, w1+j+1
                                   }, f, 9 );
    }
    //corners
    dst[0] = apply_filter( im, (int[]) {
                               0, 0, 1,
                               0, 0, 1,
                               w, w, w+1
                           }, f, 9 );
    dst[w-1] = apply_filter( im, (int[]) {
                                 w-2,   w-1,   w-1,
                                 w-2,   w-1,   w-1,
                                 2*w-2, 2*w-1, 2*w-1
                             }, f, 9 );
    dst[(h-1)*w] = apply_filter( im, (int[]) {
                                     (h-2)*w, (h-2)*w, (h-2)*w+1,
                                     (h-1)*w, (h-1)*w, (h-1)*w+1,
                                     (h-1)*w, (h-1)*w, (h-1)*w+1
                                 }, f, 9 );
    int w0 = (h-1)*w-2, w1 = (h-1)*w-1;
    dst[h*w-1] = apply_filter( im, (int[]) {
                                   w0,    w1,    w1,
                                   h*w-2, h*w-1, h*w-1,
                                   h*w-2, h*w-1, h*w-1
                               }, f, 9 );
}

/**
 * @brief Implementation of the MATLAB rgb2gray function
 *
 * See: http://www.mathworks.com/help/images/ref/rgb2gray.html
 *
 * @param gray (out) Output image
 * @param rgb Input image
 * @param npixels Size of image in pixels
 */
void rgb2gray( double* dst, double *im, size_t npixels ) {
    for(int i = 0; i < npixels; i++){
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        dst[i] = 0.2989 * r + 0.5870 * g + 0.1140 * b; //retain luminance, discard hue and saturation
    }
}


