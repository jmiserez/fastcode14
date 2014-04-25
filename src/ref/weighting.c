#include<stddef.h>
#include<math.h>
#include"weighting.h"
#include"image.h"
#include"matrix.h"
#include"convolution.h"

void weight_matrix( double* weight_matrix, double* weights, size_t len,
                    double exp ) {
    int i;
    for( i = 0; i < len; i++ )
        weight_matrix[i] *= pow(weights[i], exp);
}

void saturation( double *C, int w, int h, double *im, segments_t* mem ) {
    size_t npixels = w * h, i;
    matrix_set( C, npixels, 0.0 );

    for( i = 0; i < npixels; i++ ) {
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        double mu = (r + g + b) / 3.0;
        C[i] = sqrt(pow(r-mu,2) + pow(g-mu,2) + pow(b-mu,2)/3.0);
    }
}

void well_exposedness( double *C, int w, int h, double *im, segments_t* mem ) {
    size_t npixels = w * h, i;
    matrix_set( C, npixels, 0.0 );

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
    size_t i;
    for( i = 0; i < len; i++ ){
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

void compute_weight_matrices( double** imgs, int w, int h, int N,
                              double contrast_parm, double sat_parm,
                              double wexp_parm, segments_t* mem ) {

    size_t npixels = w * h, pix_i;
    int i;
    double* temp_matrix = mem->temp_matrix;
    double** weight_matrices = mem->weight_matrices;

    for( i = 0; i < N; i++ ) {
        if( contrast_parm > 0.0 ) {
            contrast( temp_matrix, w, h, imgs[i], mem );
            weight_matrix( weight_matrices[i], temp_matrix, npixels,
                           contrast_parm );
        }

        if( sat_parm > 0.0 ) {
            saturation( temp_matrix, w, h, imgs[i], mem );
            weight_matrix( weight_matrices[i], temp_matrix, npixels,
                           sat_parm );
        }

        if( wexp_parm > 0.0 ) {
            well_exposedness( temp_matrix, w, h, imgs[i], mem );
            weight_matrix( weight_matrices[i], temp_matrix, npixels,
                           wexp_parm );
        }

    }

    // normalize
    double sum;
    for( pix_i = 0; pix_i < npixels; pix_i++ ) {
        sum = 0.0;
        for( i = 0; i < N; i++ )
            sum += weight_matrices[i][pix_i] + 1.0E-12;
        for( i = 0; i < N; i++ )
            weight_matrices[i][pix_i] /= sum;
    }
 }
