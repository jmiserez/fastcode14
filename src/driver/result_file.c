#include "result_file.h"

FILE* res_file_open( char* path ) {
    return fopen( path, "a" );
}


void res_file_close( FILE* f ) {
    fclose(f);
}

int write_result( FILE* res_file, testconfig_t* tc, result_t* result ) {
    // TODO
    return 0;
}


