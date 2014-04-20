#ifndef IMAGE_IO_H
#define IMAGE_IO_H

double* load_tiff_rgb( uint32_t* r_width, uint32_t* r_height, const char* path );
void store_tiff_rgb( double* rgb_image, uint32_t width, uin32_t height, const char* path );
int debug_tiff_test( const char *in_img, char *out_img );

#endif // IMAGE_IO_H
