#include <fusion.h>

#include <stdio.h>
#include <string.h>
#include <stdint.h>
//#include <stdint-gcc.h> Not available on RedHat6@CAB. Not needed.
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>

#include "debug.h"

#include "testconfig.h"
#include "image_io.h"

#include "perfmon_wrapper.h"
#include "cost_model.h"

#define TIFF_DEBUG_IN "gradient.tif"
#define TIFF_DEBUG_OUT "out.tif"
#define TIFF_DEBUG_OUT2 "gradient.o.tif"

COST_VARIABLES_HERE


typedef struct {
    char* log_file;
    char* out_file; ///< output image written, if desired
    char* val_file; ///< comparision image, if desired
    double val_threshold;
} cli_options_t;

/*
 * GENERAL I/O
 **/
const char* usage_str = "./driver [options] <heights> <widths> "
        "<contrastParm> <saturationParm> <wellexpParm>\n\n"
        "options:\n"
        " --testlibtiff:\n"
        "  performs some tests using libtiff\n"
        "\n"
        " --validate <ref>\n"
        "  compare calculated result against <ref>.\n"
        "  The driver can be configured to execute exposure fusion for\n"
        "  different widths and heights. If the current width and height are\n"
        "  smaller than the width and height of the reference image, the\n"
        "  current result image is validated against the rectangular\n"
        "  subsection with respective dimensions.\n"
        "\n"
        " --threshold <double>\n"
        "  threshold for validation.\n"
        "\n"
        " --store <prefix>\n"
        "  store result into <prefix>-<w>x<h>.tif\n"
        "\n\n"
        "Example:\n"
        " $ ./driver --val val.tif --store out.tif "
        "100:50:1000\\\n100:50:1000 0.1 0.2 0.3 input1.tif input2.tif\n\n"
        "When executed using these arguments, the driver reads files \n"
        "input1.tif and input 2.tif and executes the exposure fusion \n"
        "algorithm for widths 100, 150, ..., 1000 and heights 100, 150, ...,\n"
        "1000 cropping the input images to the respective dimensions.\n\n"
        "The contrast parameter is 0.1, the saturation parameter 0.2 and the\n"
        "well-exposedness parameter is 0.3.\n\n"
        "The resulting image is written to out.tif and validated against the\n"
        "image stored in val.tif.\n";

void print_usage() {
    printf("%s", usage_str);
}

/*
 * PERFORMANCE MEASUREMENT
 **/
int run_testconfiguration( cli_options_t* cli_opts, testconfig_t* tc ) {

    int ret;

    // Initialize performance counters
    char *events[] = {
        "PERF_COUNT_HW_CPU_CYCLES",
        "CACHE-REFERENCES",
        "CACHE-MISSES",
        NULL
    };

    struct perf_data *perf_data;
    
    ret = perf_init(events, &perf_data);
    if (ret < 0) {
        FUSION_ERR("Could not initialize perfmon.\n");
        perf_cleanup(perf_data);
        return -1;
    }

    // Initialize image data
    size_t w_orig, h_orig, w, h, img_count;
    double** input_images = tc_read_input_images( &img_count, &w_orig, &h_orig,
                                                  tc );
#ifdef DEBUG
    printf("img_count: %ld, w_orig: %ld, h_orig: %ld\n", img_count, w_orig,
           h_orig);
#endif

    size_t val_read_w = 0;
    size_t val_read_h = 0;
    double *val = NULL;
    if(cli_opts->val_file != NULL){
        //load validation image
        val = load_tiff_rgb( (uint32_t*)&val_read_w,
                                   (uint32_t*)&val_read_h,
                                   cli_opts->val_file );
        assert(val != NULL);
        assert(val_read_w == w_orig);
        assert(val_read_h == h_orig);
    }

    // Perform Computations
    int err = 0;
    for( w = tc->width_range.start; w <= tc->width_range.stop && !err;
         w += tc->width_range.step ) {
        for( h = tc->height_range.start; h <= tc->height_range.stop && !err;
             h += tc->height_range.step ) {
#ifndef NDEBUG
            printf( "fusion for w: %ld, h: %ld\n", w, h );
#endif
            // crop the images
            double** target_images = crop_topleft_rgbs( input_images, w_orig,
                                                        h_orig, img_count,
                                                        w, h, true );

            void *fusion_data; ///< memory segments used by the current
            // implementation
            ret = fusion_alloc(&fusion_data, w, h, img_count);
            if (ret < 0) {
                FUSION_ERR("Error in fusion_alloc(*,w=%zu,h=%zu,N=%zu)\n",
                        w, h, img_count);
                fusion_free(fusion_data);
                err = -1;
            } else { // fusion alloc succeeded

                //TODO: warmup caches

                // Reset Counters
                COST_MODEL_RESET;
                perf_reset(perf_data);

                // Run exposure fusion
                perf_start(perf_data);
                double *result = fusion_compute(target_images, w, h, img_count, tc->contrast,
                        tc->saturation, tc->well_exposed, fusion_data);
                perf_stop(perf_data);

                // Print results
                perf_update_values(perf_data);

                //store image
                if(cli_opts->out_file != NULL){
                    char sbuf[255];
                    snprintf(sbuf,255,"%s-%dx%d.tif",cli_opts->out_file,w,h);
                    store_tiff_rgb(result,w,h,sbuf);
                }
                printf("--- Showing results for: w=%zu, h=%zu, N=%zu\n",
                        w, h, img_count);
                //validate image
                if(cli_opts->val_file != NULL){
                    //crop validation image
                    assert(val_read_w >= w);
                    assert(val_read_h >= h);
                    double *val_crop = crop_topleft_rgb(val,val_read_w,val_read_h,w,h,true);
                    assert(val_crop != NULL);

                    printf("Performance Counters:\n");

                    double rmse = rmse(result, val_crop, w, h);
                    printf("  RMSE          : %.0lf\n", rmse);
                    if(rmse > cli_opts->val_threshold){
                        FUSION_ERR("Error in validation(*,w=%zu,h=%zu,rmse=%lf)\n",
                                w, h, rmse);
                        err = -2;
                    }
                    free(val_crop);
                }

                printf("Performance Counters:\n");
                for (struct perf_data *iter = perf_data; iter->name; ++iter) {
                    printf("  %-50s : %"PRId64"\n", iter->name, iter->value);
                }

                printf("Cost Model:\n");
                printf("  Add: %"PRI_COST"\n", COST_ADD);
                printf("  Mul: %"PRI_COST"\n", COST_MUL);
                printf("  Div: %"PRI_COST"\n", COST_DIV);
                printf("  Pow: %"PRI_COST"\n", COST_POW);
                printf("  Abs: %"PRI_COST"\n", COST_ABS);

                printf("Derived Results\n");
                double flops      = COST_ADD + COST_MUL + COST_DIV;
                double cycles     = perf_data[0].value; // index as in 'events'
                double cache_load = perf_data[1].value;
                double cache_miss = perf_data[2].value;
                printf("  Flops           : %.0lf\n", flops);
                printf("  Performance     : %.3lf\n", flops/cycles);
                printf("  Cache Miss Rate : %.3lf\n", cache_miss/cache_load);

                // Cleanup
                fusion_free(fusion_data);
            } // fusion alloc succeeded

            // free target images
            free_rgbs( target_images, img_count );
        } // heights loop
    } // widths loop

    if(val != NULL){
        free(val);
    }
    // Cleanup
    tc_free_input_images( input_images, img_count );
    perf_cleanup(perf_data);

    return err;
}

int parse_cli(cli_options_t* cli_opts, testconfig_t* testconfig,
                          int argc, char* argv[]) {

    // getopt for command-line parsing. See the getopt(3) manpage
    int c;
    while (true) {
        static struct option long_options[] = {
            {"testlibtiff", no_argument,       0, 'd'},
            {"store",       required_argument, 0, 's'},
            {"validate",    required_argument, 0, 'v'},
            {"threshold",   required_argument, 0, 't'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "t", long_options, &option_index);
        if (c == -1) { // -1 indicates end of options reached
            break;
        }
        switch (c) {
        case 0:
            // the long option with name long_options[option_index].name is
            // found
            printf("getopt error on long option %s\n",
                   long_options[option_index].name);
            break;

        case 'd':
            printf("getopt: testlibtiff\n");
            debug_tiff_test( TIFF_DEBUG_IN, TIFF_DEBUG_OUT );

            uint32_t dbg_w, dbg_h;
            double *debug_rgb_image = load_tiff_rgb( &dbg_w, &dbg_h,
                                                     TIFF_DEBUG_IN );
            printf( "dbg_w: %d, dbg_h: %d\n", dbg_w, dbg_h );
            store_tiff_rgb( debug_rgb_image, dbg_w, dbg_h, TIFF_DEBUG_OUT2 );

            uint32_t* raster = rgb2tiff(debug_rgb_image, dbg_w*dbg_h);
            double err = compare_tif(raster, dbg_w, dbg_h, TIFF_DEBUG_IN );
            printf("error is: %lf\n", err);
            free_tiff( raster );
            free_rgb( debug_rgb_image );
            break;
        case 's':
            cli_opts->out_file = optarg;
            break;
        case 'v':
            cli_opts->val_file = optarg;
            break;
        case 't':
            cli_opts->val_threshold = atof(optarg);
            break;
        case '?':
            printf("getopt: error on character %c\n", optopt);
            break;

        default:
            printf("getopt: general error\n");
            abort();
        }
    }

    int num_args_remaining = argc - optind;
    int ret_optind, tmp_optind = optind;
    if( (ret_optind = read_testconfiguration(testconfig, num_args_remaining,
                                             &argv[optind])) < 0 ) {
        int err_arg = tmp_optind + abs(ret_optind) - 1;
        fprintf(stderr, "error while parsing argument %d: %s\n", err_arg+1,
               argv[err_arg]);
        return -err_arg;
    } else if ( ret_optind < 5 ) {
        fprintf( stderr, "Error: too few arguments\n" );
        return -argc;
    }
    return ret_optind + tmp_optind;
    // next time, we use a language that supports exceptions for the driver code
    // ;-)
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
    testconfig_t testconfig = {
        .width_range  = {0,0,0},
        .height_range = {0,0,0},
        .contrast     = 0.0,
        .saturation   = 0.0,
        .well_exposed = 0.0,
        .n_inputfiles = 0,
        .input_paths  = NULL,
        .ref_path     = NULL
    };

    cli_options_t cli_opts = {
        .out_file = NULL,
        .val_file = NULL
    };

    if ( parse_cli( &cli_opts, &testconfig, argc, argv ) > 0 ) {
#ifndef NDEBUG
        print_testconfiguration( &testconfig );
        if( cli_opts.out_file )
            printf("out file: %s\n", cli_opts.out_file );
        if( cli_opts.val_file )
            printf("val file: %s\n", cli_opts.val_file );
#endif
        return run_testconfiguration( &cli_opts, &testconfig );
    } else print_usage();

    return -1;
}
