#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdint-gcc.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>

#include <tiffio.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int main(int argc, char *argv[]);
void run(uint32_t **images, uint32_t nimages, uint32_t width, uint32_t height);

void exposure_fusion(double** I, int r, int c, int N, double m[3], double* R);
void contrast(double *im, uint32_t r, uint32_t c, double *C);
void saturation(double *im, uint32_t npixels, double *C);
void well_exposedness(double *im, uint32_t npixels, double *C);
void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double **pyr);

uint32_t compute_nlev(uint32_t r, uint32_t c);
void malloc_foreach(double **dst, size_t size, uint32_t N);
void malloc_gaussian_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr);
double sum3(double r, double g, double b);
double stdev3(double r, double g, double b);
void unzip(double *src, size_t src_len, uint32_t n, double **dst);
void zip(double **src, size_t src_len, uint32_t n, double *dst);
void fold(double *src, size_t src_len, double (*func)(double,double), double *dst);
void fold3(double *src, size_t src_len, double (*func)(double,double,double), double *dst);

void rgb2gray(double *im, size_t npixels, double* dst);
void ones(double *dst, size_t len);
void ones_foreach(double **dst, size_t len, uint32_t N);
void zeros(double *dst, size_t len);
void zeros_foreach(double **dst, size_t len, uint32_t N);
void scalar_add(double *src, size_t src_len, double val, double *dst);
void scalar_mult(double *src, size_t src_len, double val, double *dst);
void scalar_div(double *src, size_t src_len, double val, double *dst);
void scalar_pow(double *src, size_t src_len, double val, double *dst);
void elementwise_mult(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_div(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_add(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_sub(double *src1, size_t src1_len, double *src2, double *dst);
void elementwise_sqrt(double *src, size_t src_len, double *dst);
void elementwise_copy(double *src, size_t src_len, double *dst);
void scalar_abs(double *src, size_t src_len, double *dst);
void conv3x3_mono_replicate(double* mono, uint32_t r, uint32_t c, double* h, double* dst);

void load_images(char **path, int nimages, uint32_t **ret_stack, uint32_t *ret_widths, uint32_t *ret_heights);
void store_image(char* path, double *R, uint32_t height, uint32_t width);

void tiff2rgb(uint32_t *tiff, size_t npixels, double* ret_rgb);
int debug_tiff_test(const char *in_img, char *out_img);


int main(int argc, char *argv[]){

    //getopt for command-line parsing. See the getopt(3) manpage
    int c;
    while(true){
        static struct option long_options[] = {
            {"testlibtiff", no_argument, 0, 't'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "t", long_options, &option_index);
        if(c == -1){ // -1 indicates end of options reached
            break;
        }
        switch(c){
            case 0: // the long option with name long_options[option_index].name is found
                printf("getopt error on long option %s\n", long_options[option_index].name);
                break;
            case 't':
                printf("getopt: testlibtiff\n");
                debug_tiff_test("gradient.tif", "out.tif");
                break;
            case '?':
                printf("getopt: error on character %c\n", optopt);
                break;
            default:
                printf("getopt: general error\n");
                abort();
        }
    }
    int num_opts = optind-1;
    int num_args_remaining = argc-optind;

    if(num_opts == 0 && num_args_remaining == 0){
        printf("Usage: ./fusion <options> <paths of images>\n");
        return 0;
    }
    if(num_args_remaining > 0){ //get rest of arguments (optind is defined in getopt.h and used by getopt)
        //use arguments

        //load all images specified on the command line
        //TODO: extract to function
        uint32_t nimages = num_args_remaining;

        assert(nimages > 0);

        char** argv_start = &argv[optind];
        uint32_t **images = malloc(nimages*sizeof(uint32_t*));
        uint32_t *image_widths = malloc(nimages*sizeof(uint32_t));
        uint32_t *image_heights = malloc(nimages*sizeof(uint32_t));
        load_images(argv_start, nimages, images, image_widths, image_heights);

        assert(images != NULL);

#ifndef NDEBUG
        for(int i = 0; i < nimages; i++){
            assert(image_widths[i] == image_widths[0]);
            assert(image_heights[i] == image_heights[0]);
        }
#endif

        run(images, nimages, image_widths[0], image_heights[0]);

        for(int i = 0; i < nimages; i++){
            free(images[i]);
        }
        free(images);
        free(image_widths);
        free(image_heights);
    }
    return 0;
}

/**
 * @brief Run validation and performance benchmarks
 */
void run(uint32_t **images, uint32_t nimages, uint32_t width, uint32_t height){
    // convert raw images to something we can work with
    uint32_t npixels = width*height;
    //malloc space for the array of pointers to the converted images
    double **I = malloc(nimages*sizeof(double*));
    for(int i = 0; i < nimages; i++){
        //malloc space for the double image
        double *converted_image = malloc(3*npixels*sizeof(double));
        assert(converted_image != NULL);

        tiff2rgb(images[i], npixels, converted_image);
        I[i] = converted_image;
        scalar_mult(I[i],npixels*3,1.0/255.0,I[i]);
    }

    // malloc space for the fused image
    size_t R_len = npixels*3;
    double *R = malloc(R_len*sizeof(double));
    assert(R != NULL);

    //TODO: make these parameters changeable

    double m[3];
    m[0] = 0.5;
    m[1] = 0.5;
    m[2] = 0.5;

    //run fusion
    exposure_fusion(I, height, width, nimages, m, R);

//    store_image("debug.tif", R, height, width);

    free(R);
    for(int i = 0; i < nimages; i++){
        free(I[i]);
    }
    free(I);
}

//
// Exposure Fusion functionality
//

/**
 * @brief Implementation of exposure_fusion.m
 * @param I represents a stack of N color images (at double precision).
 *        Dimensions are (height x width x 3 x N).
 * @param r height of image
 * @param c width of image
 * @param N number of images
 * @param m 3-tuple that controls the per-pixel measures. The elements
 *        control contrast, saturation and well-exposedness,
 *        respectively.
 */
void exposure_fusion(double** I, int r, int c, int N, double m[3], double* R){
    size_t I_len = N;
    size_t I_len2 = r*c*3;
    size_t npixels = r*c;

    double contrast_parm = m[0];
    double sat_parm = m[1];
    double wexp_parm = m[2];

    //W[n] is a weight map (1 value/pixel)
    //There is one for each of the N images
    size_t W_len = N;
    size_t W_len2 = npixels;
    double** W = malloc(W_len*sizeof(double*));
    assert(W != NULL);
    assert(W_len == I_len);

    for (int n = 0; n < N; n++){
        W[n] = malloc(W_len2*sizeof(double));
        assert(W[n] != NULL);
    }

    //C is used as a temporary variable (1 value/pixel)
    size_t C_len = npixels;
    double* C = malloc(C_len*sizeof(double));
    assert(C != NULL);
    assert(W_len2 == C_len);

    //for each image, calculate the weight maps
    for (int n = 0; n < N; n++){
        ones(W[n], W_len2);

        if(contrast_parm > 0){
            contrast(I[n],r,c,C);
            scalar_pow(C,C_len,contrast_parm,C);
            elementwise_mult(W[n],W_len2,C,W[n]);
        }

        if(sat_parm > 0){
            saturation(I[n],npixels,C);
            scalar_pow(C,C_len,sat_parm,C);
            elementwise_mult(W[n],W_len2,C,W[n]);
        }

        if(wexp_parm > 0){
            well_exposedness(I[n],npixels,C);
            scalar_pow(C,C_len,wexp_parm,C);
            elementwise_mult(W[n],W_len2,C,W[n]);
        }

        scalar_add(W[n],W_len2,1.0E-12,W[n]);

    }

    //normalize weights: the total sum of weights for each pixel should be 1 across all N images
    elementwise_copy(W[0],W_len2,C);
    for (int n = 1; n < N; n++){
        elementwise_add(C,C_len,W[n],C);
    }
    for (int n = 0; n < N; n++){
        elementwise_div(W[n],W_len2,C,W[n]);
    }
#ifndef NDEBUG
    for(int i = 0; i < W_len2; i++){
        double sum_weight = 0;
        for (int n = 0; n < W_len; n++){
            sum_weight += W[n][i];
        }
        assert(sum_weight >= 0.99 && sum_weight <= 1.01); //ensure all weights sum to one for each pixel
    }
#endif
    uint32_t nlev = compute_nlev(r,c);
    assert(nlev != 0);

    //create empty pyramid
    double **pyr = NULL;
    malloc_gaussian_pyramid(r,c,3,nlev,&pyr);
    assert(pyr != NULL);



    //TODO

    free(C);
    for(int i = 0; i < nlev; i++){
        free(pyr[i]);
    }

    free(pyr);

    //TODO: remove this and replace with fused image
    //output image = 2nd input image
    elementwise_copy(I[0],I_len2,R);

    printf("done\n");

    for (int n = 0; n < N; n++){
        free(W[n]);
    }
    free(W);
}

void contrast(double *im, uint32_t r, uint32_t c, double *C){
    //laplacian filter
    double h[] = {
        0.0, 1.0, 0.0,
        1.0, -4.0, 1.0,
        0.0, 1.0, 0.0
    };
    zeros(C, r*c);
    //for each image, calculate contrast measure on grayscale version of the image
    size_t mono_len = r*c; //1 value/pixel
    double *mono = malloc(mono_len*sizeof(double));
    assert(mono != NULL);
    rgb2gray(im, r*c, mono);
    conv3x3_mono_replicate(mono,r,c,h,C);
    free(mono);
}


void saturation(double *im, uint32_t npixels, double *C){
    //saturation is computed as the standard deviation of the color channels
    size_t C_len = npixels;
    zeros(C, C_len);

    // simple version
    for(int i = 0; i < npixels; i++){
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        double mu = (r + g + b) / 3.0;
        C[i] = sqrt(pow(r-mu,2) + pow(g-mu,2) + pow(b-mu,2)/3.0);
    }

    // Matlab-like literal version (operations on whole matrices)
    //
    //    size_t R_len = npixels; //1 value/pixel
    //    double* R = malloc(R_len*sizeof(double));
    //    assert(R != NULL);
    //    size_t G_len = npixels; //1 value/pixel
    //    double* G = malloc(G_len*sizeof(double));
    //    assert(G != NULL);
    //    size_t B_len = npixels; //1 value/pixel
    //    double* B = malloc(B_len*sizeof(double));
    //    assert(B != NULL);

    //    double* rgb[] = {
    //        R, G, B
    //    };

    //    unzip(im,3*npixels,3,rgb); //split image into 3 channels

    //    size_t mu_len = npixels; //1 value/pixel
    //    double* mu = malloc(R_len*sizeof(double));
    //    assert(mu != NULL);

    //    elementwise_add(R,R_len,G,mu);
    //    elementwise_add(mu,mu_len,B,mu);
    //    scalar_div(mu,mu_len,3,mu);

    //    elementwise_sub(R,R_len,mu,R);
    //    elementwise_sub(G,G_len,mu,G);
    //    elementwise_sub(B,B_len,mu,B);

    //    scalar_pow(R,R_len,2,R);
    //    scalar_pow(G,G_len,2,G);
    //    scalar_pow(B,B_len,2,B);

    //    elementwise_add(R,R_len,G,C);
    //    elementwise_add(C,C_len,B,C);
    //    scalar_div(C,C_len,3,C);

    //    free(R);
    //    free(G);
    //    free(B);
    //    free(mu);
}

void well_exposedness(double *im, uint32_t npixels, double *C){
    size_t C_len = npixels;
    zeros(C, C_len);

    double sig = 0.2;
    for(int i = 0; i < npixels; i++){
        double r = im[i*3];
        double g = im[i*3+1];
        double b = im[i*3+2];
        r = exp(-0.5*pow(r - 0.5,2) / pow(sig,2));
        g = exp(-0.5*pow(g - 0.5,2) / pow(sig,2));
        b = exp(-0.5*pow(b - 0.5,2) / pow(sig,2));
        C[i] = r*g*b;
    }
}

void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t nlev, double **pyr){
    //pyr is an array of nlev arrays containing images of different (!) sizes
    //at this point pyr is already malloc-ed

    //TODO

}

//
// Helper functions
//

/**
 * Compute the highest possible pyramid
 */
uint32_t compute_nlev(uint32_t r, uint32_t c){
    return (uint32_t)(floor((log2(MIN(r,c)))));
}

void malloc_foreach(double **dst, size_t size, uint32_t N){
    for(int i = 0; i < N; i++){
        dst[i] = malloc(size);
        assert(dst[i] != NULL);
    }
}

/**
 * Allocate memory for the gaussian pyramid at *pyr
 */
void malloc_gaussian_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr){
    size_t pyr_len = nlev;
    *pyr = (double**) malloc(pyr_len*sizeof(double*));
    assert(*pyr != NULL);

    uint32_t r_level = r;
    uint32_t c_level = c;

    for(int i = 0; i < nlev; i++){
        // width if odd: (W-1)/2+1, otherwise: (W-1)/2
        r_level = (r_level % 2 == 0) ? (r_level-1)/2 : (r_level-1)/2+1;
        c_level = (c_level % 2 == 0) ? (c_level-1)/2 : (c_level-1)/2+1;

        size_t L_len = r_level*c_level*channels;
        double* L = malloc(L_len*sizeof(double));
        assert(L != NULL);

        (*pyr)[i] = L; //add entry to array of pointers to image levels

    }
}

//
// The following functions are for debugging only.
//

double sum3(double r, double g, double b){
    return r+g+b;
}
double stdev3(double r, double g, double b){
    double mu = (r+g+b)/3.0;
    return sqrt((pow(r-mu,2) + pow(g-mu,2) + pow(b-mu,2))/3.0);
}

/**
 * @brief Unzip an array of n-tuples into n arrays of 1 tuples
 */
void unzip(double *src, size_t src_len, uint32_t n, double **dst){
    assert(src_len % n == 0); //ensure there are no left-over elements
    for(int i = 0; i < src_len/n; i++){
        for(int k = 0; k < n; k++){
            dst[k][i] = src[i*n+k];
        }
    }
}

/**
 * @brief Zip n arrays of 1-tuples into 1 array of n-tuples
 */
void zip(double **src, size_t src_len, uint32_t n, double *dst){
    for(int i = 0; i < src_len; i++){
        for(int k = 0; k < n; k++){
            dst[i*n+k] = src[k][i];
        }
    }
}

/**
 * @brief Sum an array
 */
void sum(double *src, size_t src_len, double *dst){
    double out = 0.0;
    for(int i = 0; i < src_len; i++){
        out += src[i];
    }
    *dst = out;
}

/**
 * @brief Folds an array of n 3-tuples into an array of n 1-tuples
 */
void fold3(double *src, size_t src_len, double (*func)(double,double,double), double *dst){
    assert(src_len % 3 == 0); //ensure there are no left-over elements
    for(int i = 0; i < src_len; i++){
        dst[i] = func(src[3*i],src[3*i+1],src[3*i+2]);
    }
}

//
// MATLAB-equivalent functionality
//

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
    }
}

void ones(double *dst, size_t len){
    for(int i = 0; i < len; i++){
        dst[i] = (double)1.0;
    }
}
void ones_foreach(double **dst, size_t len, uint32_t N){
    for(int i = 0; i < N; i++){
        ones(dst[i], len);
    }
}

void zeros(double *dst, size_t len){
    for(int i = 0; i < len; i++){
        dst[i] = (double)0.0;
    }
}
void zeros_foreach(double **dst, size_t len, uint32_t N){
    for(int i = 0; i < N; i++){
        zeros(dst[i],len);
    }
}

void scalar_add(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i] + val;
    }
}

void scalar_mult(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i] * val;
    }
}

void scalar_div(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i] / val;
    }
}

void scalar_pow(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = pow(src[i],val);
    }
}

void elementwise_mult(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] * src2[i];
    }
}

void elementwise_div(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] / src2[i]; //beware of division by zero
    }
}

void elementwise_add(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] + src2[i];
    }
}

void elementwise_copy(double *src, size_t src_len, double *dst){
//    memcpy(dst,src,src_len*sizeof(double));
    printf("%p %p\n",src,dst);
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i];
    }
}

void elementwise_sub(double *src1, size_t src1_len, double *src2, double *dst){
    for(int i = 0; i < src1_len; i++){
        dst[i] = src1[i] - src2[i];
    }
}

void elementwise_sqrt(double *src, size_t src_len, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = sqrt(src[i]);
    }
}

void scalar_abs(double *src, size_t src_len, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = fabs(src[i]);
    }
}

/**
 * @brief convolution with border replication
 */
void conv3x3_mono_replicate(double* mono, uint32_t r, uint32_t c, double* h, double* dst){
    for(int i = 1; i < r-1; i++){
        for(int j = 1; j < c-1; j++){
            dst[i*c+j] =
                    dst[(i-1)*c+(j-1)]*h[0] + dst[(i-1)*c+(j)]*h[1] + dst[(i-1)*c+(j+1)]*h[1] +
                    dst[(i)  *c+(j-1)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j+1)]*h[4] +
                    dst[(i+1)*c+(j-1)]*h[5] + dst[(i+1)*c+(j)]*h[6] + dst[(i+1)*c+(j+1)]*h[7];
        }
    }
    //edges
    for(int i = 1; i < r-1; i++){
        int j = 0;
        dst[i*c+j] =
                dst[(i-1)*c+(j)]*h[0] + dst[(i-1)*c+(j)]*h[1] + dst[(i-1)*c+(j+1)]*h[1] +
                dst[(i)  *c+(j)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j+1)]*h[4] +
                dst[(i+1)*c+(j)]*h[5] + dst[(i+1)*c+(j)]*h[6] + dst[(i+1)*c+(j+1)]*h[7];
        j = c-1;
        dst[i*c+j] =
                dst[(i-1)*c+(j-1)]*h[0] + dst[(i-1)*c+(j)]*h[1] + dst[(i-1)*c+(j)]*h[1] +
                dst[(i)  *c+(j-1)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j)]*h[4] +
                dst[(i+1)*c+(j-1)]*h[5] + dst[(i+1)*c+(j)]*h[6] + dst[(i+1)*c+(j)]*h[7];
    }
    for(int j = 1; j < c-1; j++){
        int i = 0;
        dst[i*c+j] =
                dst[(i)  *c+(j-1)]*h[0] + dst[(i)  *c+(j)]*h[1] + dst[(i)  *c+(j+1)]*h[1] +
                dst[(i)  *c+(j-1)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j+1)]*h[4] +
                dst[(i+1)*c+(j-1)]*h[5] + dst[(i+1)*c+(j)]*h[6] + dst[(i+1)*c+(j+1)]*h[7];
        i = r-1;
        dst[i*c+j] =
                dst[(i-1)*c+(j-1)]*h[0] + dst[(i-1)*c+(j)]*h[1] + dst[(i-1)*c+(j+1)]*h[1] +
                dst[(i)  *c+(j-1)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j+1)]*h[4] +
                dst[(i)  *c+(j-1)]*h[5] + dst[(i)  *c+(j)]*h[6] + dst[(i)  *c+(j+1)]*h[7];
    }
    //corners
    int i = 0;
    int j = 0;
    dst[i*c+j] =
            dst[(i)*c+(j)]*h[0] + dst[(i)*c+(j)]*h[1] + dst[(i)*c+(j+1)]*h[1] +
            dst[(i)  *c+(j)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j+1)]*h[4] +
            dst[(i+1)*c+(j)]*h[5] + dst[(i+1)*c+(j)]*h[6] + dst[(i+1)*c+(j+1)]*h[7];
    i = 0;
    j = c-1;
    dst[i*c+j] =
            dst[(i)  *c+(j-1)]*h[0] + dst[(i)  *c+(j)]*h[1] + dst[(i)  *c+(j)]*h[1] +
            dst[(i)  *c+(j-1)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j)]*h[4] +
            dst[(i+1)*c+(j-1)]*h[5] + dst[(i+1)*c+(j)]*h[6] + dst[(i+1)*c+(j)]*h[7];
    i = r-1;
    j = 0;
    dst[i*c+j] =
            dst[(i-1)*c+(j-1)]*h[0] + dst[(i-1)*c+(j)]*h[1] + dst[(i-1)*c+(j+1)]*h[1] +
            dst[(i)  *c+(j-1)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j+1)]*h[4] +
            dst[(i)  *c+(j-1)]*h[5] + dst[(i)  *c+(j)]*h[6] + dst[(i)  *c+(j+1)]*h[7];
    i = r-1;
    j = c-1;
    dst[i*c+j] =
            dst[(i-1)*c+(j-1)]*h[0] + dst[(i-1)*c+(j)]*h[1] + dst[(i-1)*c+(j)]*h[1] +
            dst[(i)  *c+(j-1)]*h[2] + dst[(i)  *c+(j)]*h[3] + dst[(i)  *c+(j)]*h[4] +
            dst[(i)  *c+(j-1)]*h[5] + dst[(i)  *c+(j)]*h[6] + dst[(i)  *c+(j)]*h[7];

}

//
// TIFF functionality
//

/**
 * @brief Load a series of images
 * @param path A list of image paths to load
 * @param nimages Number of images to load
 * @param reduce Factor by which to reduce the image sizes
 * @param ret_stack (out) Stack of images
 */
void load_images(char **path, int nimages, uint32_t **ret_stack, uint32_t *ret_widths, uint32_t *ret_heights){
    for(int i = 0; i < nimages; i++){
        TIFF* tif = TIFFOpen(path[i], "r");
        assert(tif != NULL);
        uint32_t width, height;
        uint32_t* raster;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
        size_t npixels = width * height;

        ret_widths[i] = width;
        ret_heights[i] = height;

        raster = (uint32_t*) _TIFFmalloc(npixels * sizeof (uint32_t));
        assert(raster != NULL);
        uint32_t *ret_img = (uint32_t*) malloc(npixels * sizeof (uint32_t));
        assert(ret_img != NULL);

        if (TIFFReadRGBAImageOriented(tif, width, height, raster, ORIENTATION_TOPLEFT, 0)) { //same as any image viewer
            for(int k = 0; k < npixels; k++){
                ret_img[k] = raster[k];
            }
        }
        ret_stack[i] = ret_img;
        _TIFFfree(raster);
        TIFFClose(tif);
    }
}

void store_image(char* path, double *R, uint32_t height, uint32_t width){
    uint32_t npixels = width*height;

    TIFF *out = TIFFOpen(path, "w");
    uint32_t* raster = (uint32_t*) _TIFFmalloc(npixels * sizeof(uint32_t));

    assert(raster != NULL);
    for(int i = 0; i < npixels; i++){
        uint32_t packed = 0;
        packed |= (uint32_t)(fmin(fmax(0.0, round(R[i*3]*255.0)), 255.0));
        packed |= ((uint32_t)(fmin(fmax(0.0, round(R[i*3+1]*255.0)), 255.0))) << 8;
        packed |= ((uint32_t)(fmin(fmax(0.0, round(R[i*3+2]*255.0)), 255.0))) << 16;
        packed |= ((uint32_t)255) << 24;
        raster[i] = packed;
    }

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, height);
    TIFFWriteEncodedStrip(out, 0, raster, width*height*sizeof(uint32_t));
    _TIFFfree(raster);
    TIFFClose(out);
}

/**
 * @brief Converts a TIFF image to an RGB array of doubles, omitting the A channel.
 * @param tiff The TIFF image in ABGR format (MSB to LSB)
 * @param pixels Size of image in pixels
 * @param rgb (out) Output image, array of npixels*3 doubles
 */
void tiff2rgb(uint32_t *tiff, size_t npixels, double* ret_rgb){
    for(int i = 0; i < npixels; i++){
        unsigned char r = TIFFGetR(tiff[i]);
        unsigned char g = TIFFGetG(tiff[i]);
        unsigned char b = TIFFGetB(tiff[i]);
        ret_rgb[i*3] = r;
        ret_rgb[i*3+1] = g;
        ret_rgb[i*3+2] = b;
    }
}
/**
 * @brief Converts a TIFF image to an RGB array of unsigned 8 bit ints, omitting the A channel.
 * @param tiff
 * @param npixels
 * @param ret_rgb
 */
void tiff2rgb8(uint32_t *tiff, size_t npixels, uint8_t* ret_rgb){
    for(int i = 0; i < npixels; i++){
        unsigned char r = TIFFGetR(tiff[i]);
        unsigned char g = TIFFGetG(tiff[i]);
        unsigned char b = TIFFGetB(tiff[i]);
        ret_rgb[i*3] = r;
        ret_rgb[i*3+1] = g;
        ret_rgb[i*3+2] = b;
    }
}
/**
 * @brief Separates a TIFF image into it's channels, omitting the A channel.
 * @param tiff
 * @param npixels
 * @param ret_r
 * @param ret_g
 * @param ret_b
 */
void tiff2channels(uint32_t *tiff, size_t npixels, double* ret_r, double* ret_g, double* ret_b){
    for(int i = 0; i < npixels; i++){
        unsigned char r = TIFFGetR(tiff[i]);
        unsigned char g = TIFFGetG(tiff[i]);
        unsigned char b = TIFFGetB(tiff[i]);
        ret_r[i] = r;
        ret_g[i] = g;
        ret_b[i] = b;
    }
}
/**
 * @brief Separates a TIFF image into it's channels, omitting the A channel.
 * @param tiff
 * @param npixels
 * @param ret_r
 * @param ret_g
 * @param ret_b
 */
void tiff2channels8(uint32_t *tiff, size_t npixels, uint8_t* ret_r, uint8_t* ret_g, uint8_t* ret_b){
    for(int i = 0; i < npixels; i++){
        unsigned char r = TIFFGetR(tiff[i]);
        unsigned char g = TIFFGetG(tiff[i]);
        unsigned char b = TIFFGetB(tiff[i]);
        ret_r[i] = r;
        ret_g[i] = g;
        ret_b[i] = b;
    }
}

int debug_tiff_test(const char *in_img, char *out_img){
    TIFF* tif = TIFFOpen(in_img, "r");

    if (tif) {
        uint32_t w, h;
        size_t npixels;
        uint32_t* raster;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
        npixels = w * h;
        printf("w=%d h=%d\n",w,h);
        raster = (uint32_t*) _TIFFmalloc(npixels * sizeof (uint32_t));
        if (raster != NULL) {
            if (TIFFReadRGBAImageOriented(tif, w, h, raster, ORIENTATION_TOPLEFT, 0)) { //same as GIMP and others
                for(int i = 0; i < npixels; i++){
                    unsigned char r = TIFFGetR(raster[i]);
                    unsigned char g = TIFFGetG(raster[i]);
                    unsigned char b = TIFFGetB(raster[i]);
                    unsigned char a = TIFFGetA(raster[i]);
                    printf("[%u %u %u %u] ", r,g,b,a);
                    if(i%w == w-1){
                        printf("\n");
                    }
                }
            }
            _TIFFfree(raster);
        }
        TIFFClose(tif);
        TIFF *out= TIFFOpen(out_img, "w");
        int samplesperpixel = 4;
        int width = 256;
        int height = 256;
        char image[width*height*samplesperpixel];
        for(int i = 0; i < 256; i++){
            for(int k = 0; k < 256*4; k+=4){
                image[i*256*4+k+0] = ((i+k)/6+k-i) % 256; //R
                image[i*256*4+k+1] = ((i+k)/4+i-k) % 256; //G
                image[i*256*4+k+2] = (i+k)/2+k % 256;  //B
                image[i*256*4+k+3] = i; //A
            }
        }
        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, samplesperpixel);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

        TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, height);
        TIFFWriteEncodedStrip(out, 0, &image[0], width*height*samplesperpixel);
        TIFFClose(out);
    }
    return 0;
}
