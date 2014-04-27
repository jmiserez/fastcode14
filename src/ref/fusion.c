#include<fusion.h>
#include<stddef.h>
#include"segments.h"
#include"matrix.h"
#include"weighting.h"
#include"image.h"
#include"convolution.h"

double* exposure_fusion(double** I, int w, int h, int N,
                        double contrast_parm, double sat_parm, double wexp_parm,
                        void* _segments) {

    size_t npixels = w * h;

    segments_t* mem = (segments_t*) _segments;
    compute_weight_matrices( I, w, h, N, contrast_parm, sat_parm, wexp_parm,
                             mem );


}



