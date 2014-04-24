#ifndef SEGMENTS_H
#define SEGMENTS_H

/*
 * The following struct contains pointers to the memory segments allocated by
 * the alloc_fusion() function.
 */
typedef struct {
    double** weight_matrices;
    double* mono_matrix;
    double* temp_matrix;
    double* result;
} segments_t;

#endif // SEGMENTS_H
