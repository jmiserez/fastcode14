#ifndef WEIGHTING_H
#define WEIGHTING_H
void compute_weight_matrices( double** imgs, int w, int h, int N,
                              double contrast_parm, double sat_parm,
                              double wexp_parm, segmets_t* mem );

#endif // WEIGHTING_H
