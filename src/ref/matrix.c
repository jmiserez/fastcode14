#include"matrix.h"

void matrix_set( double *dst, size_t len, double val ) {
    size_t i;
    for( i = 0; i < len; i++ )
        dst[i] = val;
}
