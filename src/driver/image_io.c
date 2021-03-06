/*
 * image_io.c
 * ----------
 *
 *
 */

#include"image_io.h"
#include<assert.h>
#include<stdlib.h>
#include<stdint.h>
#include<tiffio.h>
#include<math.h>

/**
 * @brief load_image loads a tiff image located at path and stores the result in
 * an array of uint32_t.
 * @param ret_width If the image can be loaded successfully, holds the width of
 * the image.
 * @param ret_height If the image can be loaded successfully, holds the height
 * of the image.
 * @param path Path of the tif image to be loaded.
 * @return A pointer to an uint32_t array holding the loaded image.
 */
uint32_t* load_image( uint32_t *ret_width, uint32_t *ret_height, char *path ) {
    TIFF* tif = TIFFOpen(path, "r");
    assert(tif != NULL);
    uint32_t *raster, *ret_value = NULL;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, ret_width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, ret_height);
    size_t npixels = (*ret_width) * (*ret_height);

    raster = (uint32_t*) _TIFFmalloc(npixels * sizeof (uint32_t));
    assert(raster != NULL);
    ret_value = raster;

    if (!TIFFReadRGBAImageOriented(tif, *ret_width, *ret_height, raster,
                                   ORIENTATION_TOPLEFT, 0)) {
        //same as any image viewer
        ret_value = NULL;
    }
    TIFFClose(tif);
    return ret_value;
}

/**
 * @brief tiff2rgb converts an image given as a rgba tif image into a double
 * array, ommitting the alpha channel.
 * @param r_rgb pointer to a double array where the resulting image shall be
 * stored.
 * @param tiff the image to be converted.
 * @param npixels the number of pixels (width*height).
 * @param norm normalization factor (will be set to default (255.0), if 0.0 is
 * given)
 */
void tiff2rgb( double *r_rgb, uint32_t *tiff, size_t npixels,
               const double _norm ){
    double norm = _norm;
    if( norm == 0.0 )
        norm = 255.0;
    for( int i = 0; i < npixels; i += 1 ){
        unsigned char r = TIFFGetR(tiff[i]);
        unsigned char g = TIFFGetG(tiff[i]);
        unsigned char b = TIFFGetB(tiff[i]);
        r_rgb[3*i  ] = r/norm;
        r_rgb[3*i+1] = g/norm;
        r_rgb[3*i+2] = b/norm;
    }
}

/**
 * @brief rgb2tiff converts a double array holding an rgb image into a tiff
 * raster, ignoring the alpha channel
 * @param rgb_image double array holding the rgb image
 * @param npixels number of pixels
 * @return raster
 */
uint32_t* rgb2tiff( const double* rgb_image, size_t npixels ) {
    uint32_t* raster = (uint32_t*) _TIFFmalloc(npixels * sizeof(uint32_t));
    for( int i = 0; i < npixels; i += 1 ) {
        uint32_t packed = 0;
        packed |= ( uint32_t)(fmin(fmax(0.0, round(rgb_image[3*i  ]*255.0)),
                              255.0));
        packed |= ((uint32_t)(fmin(fmax(0.0, round(rgb_image[3*i+1]*255.0)),
                              255.0))) << 8;
        packed |= ((uint32_t)(fmin(fmax(0.0, round(rgb_image[3*i+2]*255.0)),
                              255.0))) << 16;
        packed |= ((uint32_t)255) << 24;
        raster[i] = packed;
    }
    return raster;
}

/**
 * @brief load_tiff_rgb loads a tif image and converts the image to a double
 * array (rgb values are normalized)
 * @param r_width width of the loaded image
 * @param r_height height of the loaded image
 * @param path path to the tif file to be loaded
 * @return pointer to the double array containing the image
 */
double* load_tiff_rgb( uint32_t* r_width, uint32_t* r_height, char* path ) {
    uint32_t *tif_img;
    double *rgb_image = NULL;
    if( (tif_img = load_image(r_width, r_height, path)) != NULL ) {
        uint32_t npixels = (*r_width)*(*r_height);
        rgb_image = (double*) malloc(npixels*3*sizeof(double));
        tiff2rgb( rgb_image, tif_img, npixels, 255.0 );
        _TIFFfree(tif_img);
    }
    return rgb_image;
}

/**
 * @brief store_tiff_raster stores a image raster, expecting 8 bits per sample
 * and 4 samples per pixel
 * @param path path to the file the image shall be written to.
 * @param raster the raster.
 * @param width width of the image.
 * @param height height of the image
 */
void store_tiff_raster( char* path, uint32_t *raster, uint32_t width,
                        uint32_t height ) {
    TIFF *out = TIFFOpen(path, "w");

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
 * @brief store_tiff_rgb stores an image given as a double array (assuming
 * normalized values) to a file.
 * @param rgb_image pointer to a double array containing the image.
 * @param width width of the image.
 * @param height height of the image.
 * @param path Path where the image shall be stored.
 *
 * Remark: Some way of returning an error?
 */
void store_tiff_rgb( double* rgb_image, uint32_t width,
                     uint32_t height, char* path ) {
    size_t npixels = width*height;
    uint32_t* raster = rgb2tiff( rgb_image, npixels );
    assert(raster != NULL);
    store_tiff_raster( path, raster, width, height );
    _TIFFfree(raster);
}

/**
 * @brief debug_tiff_test
 * @param in_img
 * @param out_img
 * @return
 */
int debug_tiff_test( const char *in_img, char *out_img ){
    TIFF* tif = TIFFOpen(in_img, "r");
    uint32_t w, h;
    size_t npixels;
    uint32_t* raster;

    if (tif) {
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
        npixels = w * h;
        printf("w=%d h=%d\n",w,h);
        raster = (uint32_t*) _TIFFmalloc(npixels * sizeof (uint32_t));
        if (raster != NULL) {
            if (TIFFReadRGBAImageOriented(tif, w, h, raster,
                                          ORIENTATION_TOPLEFT, 0)) {
                //same as GIMP and others
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
        store_tiff_raster(out_img, (uint32_t*) image, width, height);
    }

    return 0;
}

/**
 * Use pixel color values from 0 to 255 for more intuitive numbers.
 * Also round values as we would when actually writing them out.
 */
double compare_rmse(double *image, double *reference, uint32_t w, uint32_t h, int *ret_num_differing_pixels){
    double res = 0.0;
#ifdef DEBUG
    printf("w: %d,h: %d\n", w, h);
#endif
    uint32_t npixels = w*h;
    int differing_pixels = 0;
    for(int i = 0; i < npixels; i++){
        double r = (fmin(fmax(0.0, round(image[i*3]*255.0)), 255.0))
                - (fmin(fmax(0.0, round(reference[i*3]*255.0)), 255.0));
        double g = (fmin(fmax(0.0, round(image[i*3+1]*255.0)), 255.0))
                - (fmin(fmax(0.0, round(reference[i*3+1]*255.0)), 255.0));
        double b = (fmin(fmax(0.0, round(image[i*3+2]*255.0)), 255.0))
                - (fmin(fmax(0.0, round(reference[i*3+2]*255.0)), 255.0));

        r = fabs(r);
        g = fabs(g);
        b = fabs(b);

        double diff = r*r+g*g+b*b;
        if(diff > 0){ //warning: double comparison, although we've rounded
            differing_pixels++;
        }
        res = res + diff;
    }
#ifdef DEBUG
    printf("Differing pixels: %d, Total Squared Error: %lf, Pixels: %d, MSE: %lf, RMSE: %lf\n",differing_pixels, res, npixels, res / ((double)npixels), sqrt(res / ((double)npixels)));
#endif
    res = sqrt(res / ((double)npixels));

    if(ret_num_differing_pixels != NULL){
        *ret_num_differing_pixels = differing_pixels;
    }
    return res;
}

void free_tiff( uint32_t* raster ) {
    _TIFFfree( raster );
}

void free_rgb( double* rgb_image ) {
    free(rgb_image);
}

void free_rgbs( double** rgb_images, size_t img_count ) {
    size_t i;
    for( i = 0; i < img_count; i++ ) {
        free(rgb_images[i]);
    }
    free( rgb_images );
}

/**
 * @brief crop_topleft_rgb extracts topleft rectangular box from rgb_image and
 * returns the resulting image as a newly allocated double array.
 * @param rgb_image
 * @param w_orig width of original image
 * @param h_orig height of original image
 * @param w width of target image
 * @param h height of target image
 * @param force_copy Indicates whether the image shall be copied even though
 * the resulting rectangular box has the same size as the original image.
 * @return target image
 */
double* crop_topleft_rgb( double* rgb_image, size_t w_orig, size_t h_orig,
                  size_t w, size_t h, bool force_copy ) {
    size_t npixels = w * h, i, j, trgt_idx, orig_idx;
    if( (w == w_orig) && (h == h_orig) && !force_copy )
        return rgb_image;
    if( (w <= w_orig) && (h <= h_orig) ) {
        double* trgt_image = (double*) malloc(3 * npixels * sizeof(double));
        for( i = 0; i < h; i++ ) {
            for( j = 0; j < w; j++ ) {
                trgt_idx = 3*(i*w + j);
                orig_idx = 3*(i*w_orig + j);
                trgt_image[trgt_idx  ] = rgb_image[orig_idx  ]; // R
                trgt_image[trgt_idx+1] = rgb_image[orig_idx+1]; // G
                trgt_image[trgt_idx+2] = rgb_image[orig_idx+2]; // B
            }
        }
        return trgt_image;
    }
    return NULL;
}

/**
 * @brief crop_topleft_rgbs applies crop_topleft_rgb to an array of images.
 * @param rgb_images source images
 * @param w_orig original width of the images
 * @param h_orig original height of the images
 * @param image_count number of images
 * @param w target width
 * @param h target height
 * @param force_copy indicates whether the source images shall be copied even
 * though they have the same size as the target images (copy if true)
 * @return pointer to an array of pointers to the images.
 */
double** crop_topleft_rgbs( double** rgb_images, size_t w_orig,
                            size_t h_orig, size_t image_count,
                            size_t w, size_t h, bool force_copy ) {
    double** trgt_images = (double**) malloc( image_count * sizeof( double ));
    size_t i;
    for( i = 0; i < image_count; i++ ) {
        trgt_images[i] = crop_topleft_rgb( rgb_images[i], w_orig, h_orig,
                                        w, h, force_copy );
    }
    return trgt_images;
}
