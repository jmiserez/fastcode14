#include<stddef.h>
#include<stdlib.h>

#include<fusion.h>
#include"segments.h"

int alloc_fusion( void* segments, int w, int h, int N ) {
    segments_t* mem = (segments_t*) segments;
    size_t npixels = w * h;
    int i;

    // allocate memory for weight matrices
    // there is one weight-matrix per image
    mem->weight_matrices = (double**) malloc( N * sizeof(double* ));
    for( i = 0; i < N; i++ ) {
        mem->weight_matrices[i] = (double*) malloc( npixels * sizeof(double));
    }

    // for the contrast weight, a monochrome version of the image
    // is needed.
    mem->mono_matrix = (double*) malloc(npixels * sizeof(double));

    // temporary matrix, also used for weight calculation.
    mem->temp_matrix = (double*) malloc(npixels * sizeof(double));

    return 1;
}

void free_fusion( void* segments, int N ) {
    int i;
    segments_t* mem = (segments_t*) segments;

    // free weight matrices
    for( i = 0; i < N; i++ )
        free( mem->weight_matrices[i] );
    free( mem->weight_matrices );

    // free monochrome matrix
    free(mem->mono_matrix);

    // free temporary matrix, used for weight calculation
    free(mem->temp_matrix);
}
