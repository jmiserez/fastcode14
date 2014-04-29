#ifndef READ_TESTCONFIG_H
#define READ_TESTCONFIG_H

#include<stdio.h>
#include<stdint.h>
#include"urange.h"



typedef struct {
    urange_t width_range;
    urange_t height_range;
    double   contrast;
    double   saturation;
    double   well_exposed;
    size_t   n_inputfiles;
    char**   input_paths;
    char*    ref_path;
} testconfig_t;

int read_testconfiguration( testconfig_t* testconfig, int argc, char* argv[] );
void print_testconfiguration( testconfig_t* tc );

double** tc_read_input_images( size_t* read_imgs, uint32_t *w, uint32_t *h,
                               testconfig_t* tc );
void tc_free_input_images( double** images, size_t n_images );

#endif // READ_TESTCONFIG_H
