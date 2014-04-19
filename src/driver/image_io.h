#ifndef IMAGE_IO_H
#define IMAGE_IO_H

void load_tiff_rgb( double* rgb_image, uint32_t* r_width,
                    uint32_t* r_height, const char* path );

void store_tiff_rgb( doubl* rgb_image, uint32_t width,
                     uin32_t height, const char* path );

#endif // IMAGE_IO_H
