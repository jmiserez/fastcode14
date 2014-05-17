#ifndef RESULT_OUT_H
#define RESULT_OUT_H

#include"testconfig.h"
#include<stdio.h>
#include<stdint.h>

typedef struct {
    size_t cycle_count;
    double error;
} result_t;


FILE* res_file_open( char* path );
void res_file_close( FILE* f );

int write_result( FILE* res_file, testconfig_t* tc, result_t* result );

#endif // RESULT_OUT_H
