
inline double apply_filter( double *im, int* idx, double* f, size_t len ) {
    size_t i;
    double res = 0.0;
    for( i = 0; i < len; i++ ) {
        res += im[idx[i]]*f[i];
    }
    return res;
}


/**
 * @brief convolution of a monochrome image with a 3x3 filter and border mode
 * "replication"
 */
void conv3x3_monochrome_replicate( double* dst, double* im, int w, int h,
                                   double* f ) {
    for(int i = 1; i < h-1; i++) {
        for(int j = 1; j < w-1; j++){
            int iw0 = (i-1)*w, iw1 = (i)*w, iw2 = (i+1)*w;
            int j0 = j-1, j1 = j, j2 = j+1;
            dst[i*w+j] = apply_filter( im, (int[]) {
                                      (iw0+j0), (iw0+j1), (iw0+j2),
                                      (iw1+j0), (iw1+j1), (iw1+j2),
                                      (iw2+j0), (iw1+j1), (iw2+j3)
                                  }, f, 9 );
        }
    }
    //edges
    for( int i = 1; i < h-1; i++ ) {
        int j = 0;
        int iw0 = (i-1)*w, iw1 = (i)*w, iw2 = (i+1)*w;

        dst[i*w+j] = apply_filter( im, (int[]) {
                                       iw0+j, iw0+j, iw0+j+1,
                                       iw1+j, iw1+j, iw1+j+1,
                                       iw2+j, iw2+j, iw2+j+1
                                   }, f, 9);
        j = w-1;
        dst[i*w+j] = apply_filter( im, (int[]) {
                                       iw0+(j-1), iw0+j, iw0+j,
                                       iw1+(j-1), iw1+j, iw1+j,
                                       iw2+(j-1), iw2+j, iw2+j
                                   }, f, 9 );
    }
    for(int j = 1; j < w-1; j++){
        int i = 0;
        dst[i*w+j] = apply_filter( im, (int[]) {
                                       j-1,   j,   j+1,
                                       j-1,   j,   j+1,
                                       w+j-1, w+j, w+j+1
                                   }, f, 9 );
        i = h-1;
        int w0 = (h-2)*w, w1 = (h-1)*w;
        dst[i*w+j] = apply_filter( im, (int[]) {
                                       w0+(j-1), w0+j, w0+j+1,
                                       w1+(j-1), w1+j, w1+j+1,
                                       w1+(j-1), w1+j, w1+j+1
                                   }, f, 9 );
    }
    //corners
    dst[0] = apply_filter( im, (int[]) {
                               0, 0, 1,
                               0, 0, 1,
                               w, w, w+1
                           }, f, 9 );
    dst[w-1] = apply_filter( im, (int[]) {
                                 w-2,   w-1,   w-1,
                                 w-2,   w-1,   w-1,
                                 2*w-2, 2*w-1, 2*w-1
                             }, f, 9 );
    dst[(h-1)*w] = apply_filter( im, (int[]) {
                                     (h-2)*w, (h-2)*w, (h-2)*w+1,
                                     (h-1)*w, (h-1)*w, (h-1)*w+1,
                                     (h-1)*w, (h-1)*w, (h-1)*w+1
                                 }, f, 9 );
    int w0 = (h-1)*w-2, w1 = (h-1)*w-1;
    dst[h*w-1] = apply_filter( im, (int[]) {
                                   w0,    w1,    w1,
                                   h*w-2, h*w-1, h*w-1,
                                   h*w-2, h*w-1, h*w-1
                               }, f, 9 );
}
