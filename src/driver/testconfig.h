#ifndef READ_TESTCONFIG_H
#define READ_TESTCONFIG_H

#include<stdio.h>

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

#endif // READ_TESTCONFIG_H
