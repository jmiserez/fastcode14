#include"testconfig.h"
#include"image_io.h"
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#include<ctype.h>
#include<stdio.h>

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
    int i = 0;
    size_t w = 0, h = 0, new_w, new_h;

    double** images = NULL;

    if( tc->n_inputfiles > 0 ) {
        images = (double**) malloc( tc->n_inputfiles * sizeof(double*));
        images[0] = load_tiff_rgb( &w, &h, tc->input_paths[0] );
        if( (images[0] != NULL) ) {
            for( i = 1; i < tc->n_inputfiles; i++ ) {
                images[i] = load_tiff_rgb( &new_w, &new_h, tc->input_paths[i] );
                if( !(new_w == w && new_h == h) || (images[i] == NULL) ) {
                    free(images[i]);
                    i--;
                    break;
                }
                w = new_w; h = new_h;
            }
        } else
            free( images[0] );
    }

    *ret_w = w; *ret_h = h;

    if( i < tc->n_inputfiles )
        images = realloc( images, i * sizeof(double*));
    *read_imgs = i;

    return images;
}

void tc_free_input_images( double** images, size_t n_images ) {
    int i;
    for( i = 0; i < n_images; i++ ) {
        free_rgb( images[i] );
    }
    free(images);
}
