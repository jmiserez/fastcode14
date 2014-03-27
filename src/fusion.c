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


int main(int argc, char *argv[]);
void run(uint32_t **images, uint32_t nimages, uint32_t width, uint32_t height);

void rgb2gray(double *im, size_t npixels, double* dst);
void ones(double *dst, size_t len);
void ones_foreach(double **dst, size_t len, uint32_t N);
void zeros(double *dst, size_t len);
void zeros_foreach(double **dst, size_t len, uint32_t N);
void scalar_mult(double *src, size_t src_len, double val, double *dst);
void scalar_pow(double *src, size_t src_len, double val, double *dst);
void elementwise_mult(double *src1, size_t src1_len, double *src2, double *dst);
void scalar_abs(double *src, size_t src_len, double *dst);
void conv3x3_mono_replicate(double* mono, uint32_t r, uint32_t c, double* h, double* dst);

void exposure_fusion(double** I, int r, int c, int N, double m[3], double* R);
void contrast(double *im, uint32_t r, uint32_t c, double *C);

void malloc_foreach(double **dst, size_t size, uint32_t N);

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
        uint32_t *images_width = malloc(nimages*sizeof(uint32_t));
        uint32_t *images_height = malloc(nimages*sizeof(uint32_t));
        load_images(argv_start, nimages, images, images_width, images_height);

        assert(images != NULL);

#ifndef NDEBUG
        for(int i = 0; i < nimages; i++){
            assert(images_width[i] == images_width[0]);
            assert(images_height[i] == images_height[0]);
        }
#endif

        run(images, nimages, images_width[0], images_height[0]);

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
    double *R = malloc(npixels*sizeof(double));
    assert(R != NULL);

    //TODO: make these parameters changeable

    double m[3];
    m[0] = 0.5;
    m[1] = 0.5;
    m[2] = 0.5;

    //run fusion
    exposure_fusion(I, height, width, nimages, m, R);

    store_image("debug.tif", R, height, width);

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

    double contrast_parm = m[0];
    double sat_parm = m[1];
    double wexp_parm = m[2];

    for (int n = 0; n < N; n++){
        size_t W_len = r*c; //1 value/pixel
        double* W = malloc(W_len*sizeof(double));
        assert(W != NULL);
        ones(W, W_len);

        if(contrast_parm > 0){
            size_t C_len = r*c; //1 value/pixel

            double* C = malloc(C_len*sizeof(double));
            assert(C != NULL);

            contrast(I[n],r,c,C);
            scalar_pow(C,C_len,contrast_parm,C);

            assert(W_len == C_len);
            elementwise_mult(W,W_len,C,W);
        }

    }

    //TODO: remove this and replace with fused image
    for(int i = 0; i < I_len2; i++){
        R[i] = I[1][i]; //just copy
    }
    printf("done\n");
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
}

//
// Helper functions
//

void malloc_foreach(double **dst, size_t size, uint32_t N){
    for(int i = 0; i < N; i++){
        dst[i] = malloc(size);
        assert(dst[i] != NULL);
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

void scalar_mult(double *src, size_t src_len, double val, double *dst){
    for(int i = 0; i < src_len; i++){
        dst[i] = src[i] * val;
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
