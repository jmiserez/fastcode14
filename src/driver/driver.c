#include <fusion.h>

#include <stdio.h>
#include <string.h>
#include <stdint.h>
//#include <stdint-gcc.h> Not available on RedHat6@CAB. Not needed.
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>

#include "result_file.h"
#include "testconfig.h"
#include "image_io.h"

#define TIFF_DEBUG_IN "gradient.tif"
#define TIFF_DEBUG_OUT "out.tif"
#define TIFF_DEBUG_OUT2 "gradient.o.tif"

/*
 * GENERAL I/O
 **/
const char* usage_str = "./driver [options] <testConfigs> <srcImagePath> <refSolutionPath> [<result log file>]\n"
        "options:\n"
        " --testlibtiff:\n"
        "  performs some tests using libtiff"
        " --store <file>\n"
        "  stores the result in <file>\n"
        "  (otherwise only the error is calculated and the result is abandoned)\n";

void print_usage() {
    printf("%s", usage_str);
    exit(1);
}

/*
 * PERFORMANCE MEASUREMENT
 **/
void run_testconfiguration( result_t* result, testconfig_t* tc ) {
    uint32_t w, h;
    size_t img_count;
    double** input_images = tc_read_input_images( &img_count, &w, &h, tc );
#ifdef DEBUG
    printf("img_count: %ld, w: %d, h: %d\n", img_count, w, h);
#endif
    //
    // size_t npixels = w * h;
    // ret_image = (double*) malloc( h*w*sizeof(double) );
    // alloc_fusion( h, w, img_count, &segments );
    // alloc_fusion( ... )
    // exposure_fusion( input_images, r, c, N, {tc[i].contrast, tc[i].saturation, tc[i].exposure}, ret_image,

    // convert result into tiff raster
    // uint32_t* raster = rgb2tiff( ret_image, npixels );


    // compare to reference
    // result->error = compare_tif(raster, w, h, tc->ref_path );

    // runtime and deviatation from reference solution

    // free resources
    // free_tiff( raster );
    // free( ret_image );
    tc_free_input_images( input_images, img_count );
    //free_fusion( &segments );
}

/*
 * MAIN
 **/

/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    const size_t max_test_configs = 1024;
    testconfig_t tc[max_test_configs];

    result_t result;
    int tc_count;
    int i;
    FILE* f_config;
    // double *ret_image;
    // fusion_segments_t segments;

    // getopt for command-line parsing. See the getopt(3) manpage
    int c;
    while(true){
        static struct option long_options[] = {
            {"testlibtiff", no_argument, 0, 't'},
            {"store", required_argument, 0, 's'},
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
            debug_tiff_test( TIFF_DEBUG_IN, TIFF_DEBUG_OUT );

            uint32_t dbg_w, dbg_h;
            double *debug_rgb_image = load_tiff_rgb( &dbg_w, &dbg_h, TIFF_DEBUG_IN );
            printf( "dbg_w: %d, dbg_h: %d\n", dbg_w, dbg_h );
            store_tiff_rgb( debug_rgb_image, dbg_w, dbg_h, TIFF_DEBUG_OUT2 );

            uint32_t* raster = rgb2tiff(debug_rgb_image, dbg_w*dbg_h);
            double err = compare_tif(raster, dbg_w, dbg_h, TIFF_DEBUG_IN );
            printf("error is: %lf\n", err);
            free_tiff( raster );
            free_rgb( debug_rgb_image );
            break;

        case 's':

            break;

        case '?':
            printf("getopt: error on character %c\n", optopt);
            break;

        default:
            printf("getopt: general error\n");
            abort();
        }
    }

    int num_args_remaining = argc-optind;

    if( num_args_remaining > 2 ) { //get rest of arguments (optind is defined in getopt.h and used by getopt)
        char* configFilePath = argv[argc - num_args_remaining    ];
        char* srcPath        = argv[argc - num_args_remaining + 1];
        char* refPath        = argv[argc - num_args_remaining + 2];
        char* resPath        = NULL;
        FILE* resFile = NULL;
        if( num_args_remaining > 3 ) {
            resPath = argv[argc - num_args_remaining + 3];
            resFile = res_file_open( resPath );
        }

#ifdef DEBUG
        printf("config file: %s\n", configFilePath);
        printf("src path:    %s\n", srcPath);
        printf("ref path:    %s\n", refPath);
#endif
        f_config = fopen( configFilePath, "r" );
        if( f_config ) {
            tc_count = read_testconfigurations( tc, max_test_configs, f_config, srcPath, refPath );
            fclose(f_config);

            for( i = 0; i < tc_count; i++ ) {
#ifdef DEBUG
                print_testconfiguration( &tc[i] );
#endif
                run_testconfiguration( &result, &tc[i] );
                if( resFile ) {
                    write_result( resFile, &tc[i], &result );
                }
                write_result( stdout, &tc[i], &result );
            }
            tcs_free( tc, tc_count );
        }
        if( resFile )
            res_file_close( resFile );
    } else
        print_usage();
    return 0;
}
