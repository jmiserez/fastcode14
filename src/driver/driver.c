#include <fusion.h>

#include <string.h>
#include <stdint.h>
//#include <stdint-gcc.h> Not available on RedHat6@CAB. Not needed.
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>

#include <tiffio.h>


void run(uint32_t **images, uint32_t nimages, uint32_t width, uint32_t height) {

    // convert to double arrays
    // call alloc_fusion(...)
    // start measurement
    // run fusion
    exposure_fusion(I, height, width, nimages, m, R);
    // stop measurement
    // store result image
    // free resources
        // call free_fusion(...)
}

/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {

    //getopt for command-line parsing. See the getopt(3) manpage
    int c;
    while(true){
        static struct option long_options[] = {
            {"testlibtiff", no_argument, 0, 't'},
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
            debug_tiff_test("gradient.tif", "out.tif");
            break;
        case '?':
            printf("getopt: error on character %c\n", optopt);
            break;
        default:
            printf("getopt: general error\n");
            abort();
        }
    }
    int num_opts = optind-1;
    int num_args_remaining = argc-optind;

    if(num_opts == 0 && num_args_remaining == 0){
        printf("Usage: ./fusion <options> <paths of images>\n");
        return 0;
    }
    if(num_args_remaining > 0){ //get rest of arguments (optind is defined in getopt.h and used by getopt)
        //use arguments

        //load all images specified on the command line
        //TODO: extract to function
        uint32_t nimages = num_args_remaining;

        assert(nimages > 0);

        char** argv_start = &argv[optind];
        uint32_t **images = malloc(nimages*sizeof(uint32_t*));
        uint32_t *image_widths = malloc(nimages*sizeof(uint32_t));
        uint32_t *image_heights = malloc(nimages*sizeof(uint32_t));
        // load_images(argv_start, nimages, images, image_widths, image_heights);

        assert(images != NULL);

#ifndef NDEBUG
        for(int i = 0; i < nimages; i++){
            assert(image_widths[i] == image_widths[0]);
            assert(image_heights[i] == image_heights[0]);
        }
#endif

        run(images, nimages, image_widths[0], image_heights[0]);

        for(int i = 0; i < nimages; i++){
            free(images[i]);
        }
        // free resources
        free(images);
        free(image_widths);
        free(image_heights);
    }
    return 0;
}
