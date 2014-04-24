#ifndef IMAGE_IO_H
#define IMAGE_IO_H

#include<stdint.h>
#include<stddef.h>

uint32_t* rgb2tiff( const double* rgb_image, size_t npixels );
double* load_tiff_rgb( uint32_t* r_width, uint32_t* r_height, char* path );
void free_rgb( double * rgb_image );
void store_tiff_rgb( double* rgb_image, uint32_t width,
                     uint32_t height, char* path );
double compare_tif( uint32_t *raster, uint32_t w, uint32_t h, char* path );
void free_tiff( uint32_t* raster );

int debug_tiff_test( const char *in_img, char *out_img );


#endif // IMAGE_IO_H
