#include <stdio.h>
#include <string.h>
#include <stdint.h>
//#include <stdint-gcc.h> Not available on RedHat6@CAB. Not needed.
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
void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *S, size_t S_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *S, size_t S_len, double *T, size_t T_len, double *U, size_t U_len, double *V, size_t V_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c);
void downsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *S, size_t S_len, uint32_t down_r, uint32_t down_c, double *dst);
void upsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *U, size_t U_len, double *V, size_t V_len, uint32_t up_r, uint32_t up_c, double *dst);

uint32_t compute_nlev(uint32_t r, uint32_t c);
void malloc_foreach(double **dst, size_t size, uint32_t N);
void malloc_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr, uint32_t **pyr_r, uint32_t **pyr_c);
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
void conv3x3_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* fxy, double* dst);
void conv5x5separable_symmetric(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double* dst);
void conv5x5separable_replicate(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double* dst);

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
    uint32_t *pyr_r = NULL;
    uint32_t *pyr_c = NULL;
    malloc_pyramid(r,c,3,nlev,&pyr, &pyr_r, &pyr_c);
    assert(pyr != NULL);
    assert(pyr_r != NULL);
    assert(pyr_c != NULL);

    //multiresolution blending
    double ***pyrW = malloc(N*sizeof(double**));
    uint32_t **pyrW_r = malloc(N*sizeof(double*));
    uint32_t **pyrW_c = malloc(N*sizeof(double*));
    assert(pyrW != NULL);
    assert(pyrW_r != NULL);
    assert(pyrW_c != NULL);

    //scratch space for gaussian/laplacian pyramid
    //TODO: optimize these away if possible
    size_t S_len = r*c*3;
    double* S = malloc(S_len*sizeof(double));
    assert(S != NULL);
    size_t T_len = r*c*3;
    double* T = malloc(T_len*sizeof(double));
    assert(T != NULL);
    //how much scratch space is needed for upsampling?
    uint32_t largest_upsampled_r = (r % 2 == 0) ? (((r-1)/2) + 2) * 2 : (((r-1)/2+1) + 2) * 2;
    uint32_t largest_upsampled_c = (c % 2 == 0) ? (((c-1)/2) + 2) * 2 : (((c-1)/2+1) + 2) * 2;
    size_t U_len = largest_upsampled_r*largest_upsampled_c*3;
    double* U = malloc(U_len*sizeof(double));
    assert(U != NULL);

    for (int n = 0; n < N; n++){
        //construct 1-channel gaussian pyramid from weights
        malloc_pyramid(r,c,1,nlev,&(pyrW[n]), &(pyrW_r[n]), &(pyrW_c[n]));
        assert(pyrW[n] != NULL);
        assert(pyrW_r[n] != NULL);
        assert(pyrW_c[n] != NULL);
        gaussian_pyramid(W[n],r,c,1,nlev,S,S_len,pyrW[n],pyrW_r[n],pyrW_c[n]);

        //TODO: construct 3-channel laplacian pyramid from images


        //TODO: weighted blend


    }


    //TODO: store

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
    conv3x3_monochrome_replicate(mono,r,c,h,C);
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

void gaussian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *S, size_t S_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    //pyr is an array of nlev arrays containing images of different (!) sizes
    //at this point pyr is already malloc-ed
    //pyr_r and pyr_c contain the sizes for each level

    assert(r == pyr_r[0]);
    assert(c == pyr_c[0]);

    //copy image to the finest level (note: MATLAB version is 1-indexed, here we use 0-indexing)

    elementwise_copy(im,r*c*channels,pyr[0]);

    size_t pyramid_filter_len = 5;
    double pyramid_filter[] = {.0625, .25, .375, .25, .0625};

    for(int v = 1; v < nlev; v++){
        //downsample image and store into level
        printf("nlev: %d, v: %d, pyr_r[v-1]: %d, pyr_c[v-1]: %d, pyr_r[v]: %d, pyr_c[v]: %d\n",nlev,v, pyr_r[v-1],pyr_c[v-1],pyr_r[v],pyr_c[v]);
        fflush(stdout);
        downsample(pyr[v-1],pyr_r[v-1],pyr_c[v-1],channels,pyramid_filter,pyramid_filter_len,S,S_len,pyr_r[v],pyr_c[v],pyr[v]);
    }
}

void laplacian_pyramid(double *im, uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double *S, size_t S_len, double *T, size_t T_len, double *U, size_t U_len, double *V, size_t V_len, double **pyr, uint32_t *pyr_r, uint32_t *pyr_c){
    //pyr is an array of nlev arrays containing images of different (!) sizes
    //at this point pyr is already malloc-ed
    //pyr_r and pyr_c contain the sizes for each level

    assert(r == pyr_r[0]);
    assert(c == pyr_c[0]);

    //copy image to the finest level (note: MATLAB version is 1-indexed, here we use 0-indexing)

    elementwise_copy(im,r*c*channels,pyr[0]);

    size_t pyramid_filter_len = 5;
    double pyramid_filter[] = {.0625, .25, .375, .25, .0625};

    //J = image
    elementwise_copy(im,T_len,T); //TODO: optimize this copy away, can use pointer swaps

    uint32_t S_r = r;
    uint32_t S_c = c;
    uint32_t T_r = r;
    uint32_t T_c = c;
    assert(S_r*S_c*channels <= S_len);
    assert(T_r*T_c*channels <= T_len);

    for(int v = 0; v < nlev-1; v++){
        //downsample image T further, store in S
        S_r = pyr_r[v+1];
        S_c = pyr_c[v+1];
        downsample(T,T_r,T_c,channels,pyramid_filter,pyramid_filter_len,S,S_len,S_r,S_c,S);

        assert(T_r*T_c == pyr_r[v]*pyr_c[v]);
        //upsample image S, store temporarily in pyramid
        upsample(S,S_r,S_c,channels,pyramid_filter,pyramid_filter_len,U,U_len,V,V_len,pyr_r[v],pyr_c[v],pyr[v]);

        //subtract pyramid from T, store difference (T - upsampled image) in pyramid
        elementwise_sub(T,T_r*T_c*channels,pyr[v],pyr[v]);

        T_r = S_r;
        T_c = S_c;
        //continue with downsampled image remainder
        elementwise_copy(S,S_r*S_c*channels,T); //TODO: optimize this copy away, can use pointer swaps (handle first and last cases!)
    }
    //coarsest level, residual low pass image
    elementwise_copy(T,T_len,pyr[nlev-1]);
}

void downsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *S, size_t S_len, uint32_t down_r, uint32_t down_c, double *dst){
    assert(filter_len == 5);
    assert(filter != NULL);
    assert(S != NULL);
    assert(r*c*channels <= S_len);

    //low pass filter
    conv5x5separable_symmetric(im,r,c,channels,filter,filter,S);
    //decimate, using every second entry
    // [1] -> [1]
    // [1 2] -> [1]
    // [1 2 3] -> [1 3]
    // [1 2 3 4] -> [1 3]
    // width if odd: (W-1)/2+1, otherwise: (W-1)/2
    assert(down_r > 0 && down_r <= (r-1)/2+1);
    assert(down_c > 0 && down_c <= (c-1)/2+1);
    assert((down_r-1)*2 < r);
    assert((down_c-1)*2 < c);
    for(int i = 0; i < down_r; i++){
        for(int j = 0; j < down_c; j++){
            for(int k = 0; k < channels; k++){
                assert(((i*2)*c+(j*2))*channels+k < r*c*channels); //bounds checking
                assert((i*down_c+j)*channels+k < down_r*down_c*channels); //bounds checking
                dst[(i*down_c+j)*channels+k] = S[((i*2)*c+(j*2))*channels+k];
            }
        }
    }
}

void upsample(double *im, uint32_t r, uint32_t c, uint32_t channels, double *filter, size_t filter_len, double *U, size_t U_len, double *V, size_t V_len, uint32_t up_r, uint32_t up_c, double *dst){
    assert(filter_len == 5);
    assert(filter != NULL);

    uint32_t padding = 1;

    //sizes with added 1 px border and size increase of 2x
    uint32_t r_upsampled = (r+2*padding)*2;
    uint32_t c_upsampled = (c+2*padding)*2;

    assert(U_len >= r_upsampled*c_upsampled*channels);
    assert(V_len == U_len);

    zeros(U,U_len);

    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            for(int k = 0; k < channels; k++){

                // i -> U_i (padding, then scaled)
                // 0 -> 2
                // 1 -> 4
                // 2 -> 6
                // 3 -> 8

                //i : 0 to r-1 -> 2 to 2*r (e.g. r=5, padded=7, scaled=14, 10, thus [0,1|2,3|4,5|6,7|8,9|10,11|12,13]

                U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[(i*c+j)*channels+k];
            }
        }
    }
    //top row
    int i = -1;
    for(int j = 0; j < c; j++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i+1)*c+(j  ))*channels+k];
        }
    }
    //bottom row
    i = r;
    for(int j = 0; j < c; j++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i-1)*c+(j  ))*channels+k];
        }
    }
    //left edge
    int j = -1;
    for(int i = 0; i < r; i++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i  )*c+(j+1))*channels+k];
        }
    }
    //right edge
    j = c;
    for(int i = 0; i < r; i++){
        for(int k = 0; k < channels; k++){
            U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i  )*c+(j-1))*channels+k];
        }
    }
    //corners
    for(int k = 0; k < channels; k++){
        i = -1;
        j = -1;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i+1)*c+(j+1))*channels+k];
        j = c;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i+1)*c+(j-1))*channels+k];
        i = r;
        j = -1;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i-1)*c+(j+1))*channels+k];
        j = c;
        U[((2*(i+padding))*c_upsampled+(2*(j+padding)))*channels+k] = 4*im[((i-1)*c+(j-1))*channels+k];
    }

    //interpolate
    conv5x5separable_replicate(U, r_upsampled, c_upsampled, channels, filter, filter, V);

    //remove the border and copy result
    assert(up_r <= r_upsampled-4);
    assert(up_c <= c_upsampled-4);

    for(int i = 0; i < up_r; i++){
        for(int j = 0; j < up_c; j++){
            for(int k = 0; k < channels; k++){
                dst[(i*up_c+j)*channels+k] = U[((i+2)*c_upsampled+(j+2))*channels+k];
            }
        }
    }
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
 * Allocate memory for the gaussian/laplacian pyramid at *pyr
 */
void malloc_pyramid(uint32_t r, uint32_t c, uint32_t channels, uint32_t nlev, double ***pyr, uint32_t **pyr_r, uint32_t **pyr_c){
    size_t pyr_len = nlev;
    *pyr = (double**) malloc(pyr_len*sizeof(double*));
    assert(*pyr != NULL);
    *pyr_r = (uint32_t*) malloc(nlev * sizeof(uint32_t));
    assert(*pyr_r != NULL);
    *pyr_c = (uint32_t*) malloc(nlev * sizeof(uint32_t));
    assert(*pyr_c != NULL);

    uint32_t r_level = r;
    uint32_t c_level = c;

    for(int i = 0; i < nlev; i++){
        (*pyr_r)[i] = r_level; //store dimension r at level i
        (*pyr_c)[i] = c_level; //store dimension c at level i

        size_t L_len = r_level*c_level*channels;
        double* L = malloc(L_len*sizeof(double));
        assert(L != NULL);

        (*pyr)[i] = L; //add entry to array of pointers to image levels

        // for next level, width if odd: (W-1)/2+1, otherwise: (W-1)/2
        r_level = (r_level % 2 == 0) ? (r_level-1)/2 : (r_level-1)/2+1;
        c_level = (c_level % 2 == 0) ? (c_level-1)/2 : (c_level-1)/2+1;
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
        dst[i] = pow(fabs(src[i]),val);
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
 * @brief convolution of a monochrome image with a 3x3 filter and border mode "replication"
 */
void conv3x3_monochrome_replicate(double* im, uint32_t r, uint32_t c, double* f, double* dst){
    for(int i = 1; i < r-1; i++){
        for(int j = 1; j < c-1; j++){
            dst[i*c+j] =
                    im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[1] +
                    im[(i)  *c+(j-1)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j+1)]*f[4] +
                    im[(i+1)*c+(j-1)]*f[5] + im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j+1)]*f[7];
        }
    }
    //edges
    for(int i = 1; i < r-1; i++){
        int j = 0;
        dst[i*c+j] =
                im[(i-1)*c+(j)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[1] +
                im[(i)  *c+(j)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j+1)]*f[4] +
                im[(i+1)*c+(j)]*f[5] + im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j+1)]*f[7];
        j = c-1;
        dst[i*c+j] =
                im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j)]*f[1] +
                im[(i)  *c+(j-1)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j)]*f[4] +
                im[(i+1)*c+(j-1)]*f[5] + im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j)]*f[7];
    }
    for(int j = 1; j < c-1; j++){
        int i = 0;
        dst[i*c+j] =
                im[(i)  *c+(j-1)]*f[0] + im[(i)  *c+(j)]*f[1] + im[(i)  *c+(j+1)]*f[1] +
                im[(i)  *c+(j-1)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j+1)]*f[4] +
                im[(i+1)*c+(j-1)]*f[5] + im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j+1)]*f[7];
        i = r-1;
        dst[i*c+j] =
                im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[1] +
                im[(i)  *c+(j-1)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j+1)]*f[4] +
                im[(i)  *c+(j-1)]*f[5] + im[(i)  *c+(j)]*f[6] + im[(i)  *c+(j+1)]*f[7];
    }
    //corners
    int i = 0;
    int j = 0;
    dst[i*c+j] =
            im[(i)*c+(j)]*f[0] + im[(i)*c+(j)]*f[1] + im[(i)*c+(j+1)]*f[1] +
            im[(i)  *c+(j)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j+1)]*f[4] +
            im[(i+1)*c+(j)]*f[5] + im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j+1)]*f[7];
    i = 0;
    j = c-1;
    dst[i*c+j] =
            im[(i)  *c+(j-1)]*f[0] + im[(i)  *c+(j)]*f[1] + im[(i)  *c+(j)]*f[1] +
            im[(i)  *c+(j-1)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j)]*f[4] +
            im[(i+1)*c+(j-1)]*f[5] + im[(i+1)*c+(j)]*f[6] + im[(i+1)*c+(j)]*f[7];
    i = r-1;
    j = 0;
    dst[i*c+j] =
            im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j+1)]*f[1] +
            im[(i)  *c+(j-1)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j+1)]*f[4] +
            im[(i)  *c+(j-1)]*f[5] + im[(i)  *c+(j)]*f[6] + im[(i)  *c+(j+1)]*f[7];
    i = r-1;
    j = c-1;
    dst[i*c+j] =
            im[(i-1)*c+(j-1)]*f[0] + im[(i-1)*c+(j)]*f[1] + im[(i-1)*c+(j)]*f[1] +
            im[(i)  *c+(j-1)]*f[2] + im[(i)  *c+(j)]*f[3] + im[(i)  *c+(j)]*f[4] +
            im[(i)  *c+(j-1)]*f[5] + im[(i)  *c+(j)]*f[6] + im[(i)  *c+(j)]*f[7];

}

/**
 * @brief convolution of a multi-channel image with a separable 5x5 filter and border mode "symmetric"
 */
void conv5x5separable_symmetric(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double* dst){
    //r is height (vertical), c is width (horizontal)

    //horizontal filter
    for(int i = 0; i < r; i++){ //all lines
        for(int j = 2; j < c-2; j++){
            for(int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] =
                        im[((i  )*c+(j-2))*channels+k]*fx[0] +
                        im[((i  )*c+(j-1))*channels+k]*fx[1] +
                        im[((i  )*c+(j  ))*channels+k]*fx[2] +
                        im[((i  )*c+(j+1))*channels+k]*fx[3] +
                        im[((i  )*c+(j+2))*channels+k]*fx[4];
            }
        }
        //left edge
        int j = 0; // 1 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j+1))*channels+k]*fx[0] +
                    im[((i  )*c+(j  ))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
        }
        j = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-1))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
        }
        //right edge
        j = c-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+1))*channels+k]*fx[4];
        }
        j = c-1; // [ ... -2 -1 0] 0 -1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j  ))*channels+k]*fx[3] +
                    im[((i  )*c+(j-1))*channels+k]*fx[4];
        }
    }
    //vertical filter
    for(int j = 0; j < c; j++){ //all columns
        for(int i = 2; i < r-2; i++){
            for(int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] =
                        im[((i-2)*c+(j  ))*channels+k]*fy[0] +
                        im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                        im[((i  )*c+(j  ))*channels+k]*fy[2] +
                        im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                        im[((i+2)*c+(j  ))*channels+k]*fy[4];
            }
        }
        //top edge
        int i = 0; // 1 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i+1)*c+(j  ))*channels+k]*fy[0] +
                    im[((i  )*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    im[((i+2)*c+(j  ))*channels+k]*fy[4];
        }
        i = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i-1)*c+(j  ))*channels+k]*fy[0] +
                    im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    im[((i+2)*c+(j  ))*channels+k]*fy[4];
        }
        //bottom edge
        i = r-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[4];
        }
        i = r-1; // [ ... -2 -1 0] 0 -1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i  )*c+(j  ))*channels+k]*fy[3] +
                    im[((i-1)*c+(j  ))*channels+k]*fy[4];
        }
    }
}

/**
 * @brief convolution of a multi-channel image with a separable 5x5 filter and border mode "replicate"
 */
void conv5x5separable_replicate(double* im, uint32_t r, uint32_t c, uint32_t channels, double* fx, double *fy, double* dst){
    //r is height (vertical), c is width (horizontal)

    //horizontal filter
    for(int i = 0; i < r; i++){ //all lines
        for(int j = 2; j < c-2; j++){
            for(int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] =
                        im[((i  )*c+(j-2))*channels+k]*fx[0] +
                        im[((i  )*c+(j-1))*channels+k]*fx[1] +
                        im[((i  )*c+(j  ))*channels+k]*fx[2] +
                        im[((i  )*c+(j+1))*channels+k]*fx[3] +
                        im[((i  )*c+(j+2))*channels+k]*fx[4];
            }
        }
        //left edge
        int j = 0; // 0 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j  ))*channels+k]*fx[0] +
                    im[((i  )*c+(j  ))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
        }
        j = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-1))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+2))*channels+k]*fx[4];
        }
        //right edge
        j = c-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j+1))*channels+k]*fx[3] +
                    im[((i  )*c+(j+1))*channels+k]*fx[4];
        }
        j = c-1; // [ ... -2 -1 0] 0 0
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j-2))*channels+k]*fx[0] +
                    im[((i  )*c+(j-1))*channels+k]*fx[1] +
                    im[((i  )*c+(j  ))*channels+k]*fx[2] +
                    im[((i  )*c+(j  ))*channels+k]*fx[3] +
                    im[((i  )*c+(j  ))*channels+k]*fx[4];
        }
    }
    //vertical filter
    for(int j = 0; j < c; j++){ //all columns
        for(int i = 2; i < r-2; i++){
            for(int k = 0; k < channels; k++){
                dst[(i*c+j)*channels+k] =
                        im[((i-2)*c+(j  ))*channels+k]*fy[0] +
                        im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                        im[((i  )*c+(j  ))*channels+k]*fy[2] +
                        im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                        im[((i+2)*c+(j  ))*channels+k]*fy[4];
            }
        }
        //top edge
        int i = 0; // 0 0 [0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i  )*c+(j  ))*channels+k]*fy[0] +
                    im[((i  )*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    im[((i+2)*c+(j  ))*channels+k]*fy[4];
        }
        i = 1; // -1 [-1 0 1 2 ... ]
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i-1)*c+(j  ))*channels+k]*fy[0] +
                    im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    im[((i+2)*c+(j  ))*channels+k]*fy[4];
        }
        //bottom edge
        i = r-2; // [ ... -2 -1 0 1] 1
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[3] +
                    im[((i+1)*c+(j  ))*channels+k]*fy[4];
        }
        i = r-1; // [ ... -2 -1 0] 0 0
        for(int k = 0; k < channels; k++){
            dst[(i*c+j)*channels+k] =
                    im[((i-2)*c+(j  ))*channels+k]*fy[0] +
                    im[((i-1)*c+(j  ))*channels+k]*fy[1] +
                    im[((i  )*c+(j  ))*channels+k]*fy[2] +
                    im[((i  )*c+(j  ))*channels+k]*fy[3] +
                    im[((i  )*c+(j  ))*channels+k]*fy[4];
        }
    }
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
