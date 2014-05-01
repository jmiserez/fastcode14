// dummy fusion. does nothing, helps testing the driver.

#include <stdlib.h>
#include <cost_model.h>

int fusion_alloc(void** _segments, int w, int h, int N)
{
    void *data = malloc(w*h*N*3*sizeof(double));
    *_segments = data;
    return data ? 0 : -1;
}

double* fusion_compute(double** I, int w, int h, int N,
                       double contrast_parm, double sat_parm, double wexp_parm,
                       void* _segments)
{
    double *mydata = _segments;
    size_t npixels = w * h;
    for (int i=0; i<w; ++i) {
        for (int j=0; j<h; ++j) {
            for (int n=0; n<N; ++j) {
                COST_INC_ADD(2);
                COST_INC_MUL(2);
                mydata[n * npixels + w*j + i] = n * npixels + w*j + i;
            }
        }
    }

    return mydata;
}

void fusion_free( void* _segments )
{
    if (_segments) free(_segments);
}
