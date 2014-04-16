#include <fusion.h>

#include <stdio.h>
#include <string.h>
#include <stdint.h>
//#include <stdint-gcc.h> Not available on RedHat6@CAB. Not needed.
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>

// #include <tiffio.h>

/*
 * GENERAL I/O
 **/
const char* usage_str = "./driver [options] <testConfigs> <srcImagePath> <refSolutionPath>\n";

void print_usage() {
    printf("%s", usage_str);
    exit(1);
}

/*
 * PERFORMANCE MEASUREMENT
 **/
void run(uint32_t **images, uint32_t nimages, uint32_t width, uint32_t height) {

    // convert to double arrays
    // call alloc_fusion(...)
    // start measurement
    // run fusion
    // exposure_fusion(I, height, width, nimages, m, R);
    // stop measurement
    // store result image
    // free resources
        // call free_fusion(...)
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

    if( num_args_remaining > 2 ) { //get rest of arguments (optind is defined in getopt.h and used by getopt)
        const char* configFilePath = argv[argc - num_args_remaining    ];
        const char* srcPath        = argv[argc - num_args_remaining + 1];
        const char* refPath        = argv[argc - num_args_remaining + 2];

#ifdef DEBUG
        printf("config file: %s\n", configFilePath);
        printf("src path:    %s\n", srcPath);
        printf("ref path:    %s\n", refPath);
#endif
    } else
        print_usage();
    return 0;
}
