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
} testconfig_t;

int read_testconfigurations( testconfig_t* testconfigs, size_t max_configs, FILE* f );
void print_testconfiguration( testconfig_t* tc );

double** tc_read_input_images( size_t* read_imgs, uint32_t *w, uint32_t *h, testconfig_t* tc, char* srcPath );
void tc_free_input_images( double** images, size_t n_images );

void tc_free( const testconfig_t* tc, size_t n_tc );

#endif // READ_TESTCONFIG_H
