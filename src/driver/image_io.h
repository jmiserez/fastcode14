#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include<stdint.h>
#include<stddef.h>
#include<stdbool.h>

uint32_t* rgb2tiff( const double* rgb_image, size_t npixels );
double* load_tiff_rgb( uint32_t* r_width, uint32_t* r_height, char* path );

void store_tiff_rgb( double* rgb_image, uint32_t width,
                     uint32_t height, char* path );
double compare_rmse(double *image, double *reference, uint32_t w, uint32_t h);
void free_tiff( uint32_t* raster );

int debug_tiff_test( const char *in_img, char *out_img );

double* crop_topleft_rgb( double* rgb_image, size_t w_orig, size_t h_orig,
                  size_t w, size_t h, bool force_copy );

double** crop_topleft_rgbs( double** rgb_images, size_t w_orig,
                            size_t h_orig, size_t image_count,
                            size_t w, size_t h, bool force_copy );


void free_rgbs( double** rgb_images, size_t img_count );
void free_rgb( double * rgb_image );


#endif // IMAGE_IO_H
