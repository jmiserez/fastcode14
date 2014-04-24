#ifndef READ_TESTCONFIG_H
#define READ_TESTCONFIG_H

#include<stdio.h>
#include<stdint.h>

typedef struct {
    char* filename;
    char* prefix;
    char* extension;
    size_t number_of_files;
    // parameters
    double contrast;
    double saturation;
    double exposure;
    char** input_paths;
    char* ref_path;
} testconfig_t;

int read_testconfigurations( testconfig_t* testconfigs, size_t max_configs, FILE* f, char* srcPath, char* refPath );
void print_testconfiguration( testconfig_t* tc );

double** tc_read_input_images( size_t* read_imgs, uint32_t *w, uint32_t *h, testconfig_t* tc );
void tc_free_input_images( double** images, size_t n_images );

void tcs_free( testconfig_t* tc, size_t n_tc );
void tc_free( testconfig_t* tc );

#endif // READ_TESTCONFIG_H
