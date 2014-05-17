#include"testconfig.h"
#include"image_io.h"
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#include<ctype.h>
#include<stdio.h>
#include <assert.h>


int read_testconfiguration(testconfig_t* testconfig, int argc, char* argv[]) {
    int optidx = 0;
    if( argc >= 5 ) {
        if( read_urange( &testconfig->width_range, argv[optidx++] ) < 0 )
            return -optidx;
        if( read_urange( &testconfig->height_range, argv[optidx++] ) < 0 )
            return -optidx;

        testconfig->contrast = atof( argv[optidx++] );
        testconfig->saturation = atof( argv[optidx++] );
        testconfig->well_exposed = atof( argv[optidx++] );

        testconfig->n_inputfiles = (argc-optidx);
        testconfig->input_paths = &argv[optidx];
    }
    return argc;
}

void print_testconfiguration( testconfig_t* tc ) {
    size_t i;
    printf("width range:  "); urange_print( &tc->width_range ); printf("\n");
    printf("height range: "); urange_print( &tc->height_range ); printf("\n");
    printf("contrast:     %lf\n", tc->contrast);
    printf("saturation:   %lf\n", tc->saturation);
    printf("exposure:     %lf\n", tc->well_exposed);
    if( tc->ref_path != NULL )
        printf("ref. path:    %s\n", tc->ref_path);
    else
        printf("ref. path:    N/A\n");
    printf("input paths:\n");
    for( i = 0; i < tc->n_inputfiles; i++ ) {
        printf("              %s\n", tc->input_paths[i]);
    }
}

double** tc_read_input_images( size_t* read_imgs, size_t *ret_w,
                               size_t *ret_h, testconfig_t* tc ) {
    int num_read = 0;
    size_t w = 0, h = 0;

    double** images = (double**) calloc( tc->n_inputfiles, sizeof(double*));

    for(int i = 0; i < tc->n_inputfiles; i++ ) {
        size_t read_w = 0;
        size_t read_h = 0;

        images[i] = load_tiff_rgb( (uint32_t*)&read_w,
                                   (uint32_t*)&read_h,
                                   tc->input_paths[i] );
        if(w == 0 && h == 0){
            w = read_w;
            h = read_h;
        }
        if(images[i] != NULL && w == read_w && h == read_h){
            num_read++;
        }

    }


    *ret_w = w; *ret_h = h;

    *read_imgs = num_read;

    return images;
}

void tc_free_input_images( double** images, size_t n_images ) {
    int i;
    for( i = 0; i < n_images; i++ ) {
        free_rgb( images[i] );
    }
    free(images);
}
