#include<assert.h>
#include<tiffio.h>

uint32_t* load_image( uint32_t *ret_width, uint32_t *ret_height, char *path ) {
    TIFF* tif = TIFFOpen(path, "r");
    assert(tif != NULL);
    uint32_t *raster, *ret_value = NULL;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &ret_width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &ret_height);
    size_t npixels = width * height;

    raster = (uint32_t*) _TIFFmalloc(npixels * sizeof (uint32_t));
    assert(raster != NULL);

    ret_value = raster;
    if (!TIFFReadRGBAImageOriented(tif, width, height, raster, ORIENTATION_TOPLEFT, 0)) {
        //same as any image viewer
        ret_value = NULL;
    }
    _TIFFfree(raster);
    TIFFClose(tif);
    return ret_value;
}

void load_tiff_rgb( double* rgb_image, uint32_t* r_width,
                    uint32_t* r_height, const char* path ) {
}

void store_tiff_rgb( doubl* rgb_image, uint32_t width,
                     uin32_t height, const char* path ) {

}
