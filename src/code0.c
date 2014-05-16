#ifndef NB
#error NB is not defined
#endif

double *A, *B, *C;

//The following code implements C = A B + C 
//with A, B and C being NB x NB matrices
//This is done by definition
void compute() {
    int i, j, k;
    for(i = 0; i < NB; ++i) { //NB times
      for(j = 0; j < NB; ++j) { //NB^2 times

        // Note: Optimization blocker: Memory aliasing. Thus I use the "out"
        //       variable, same as in the main.c example. The direct approach
        //       C[i*NB+j] += A[i*NB+k] * B[k*NB+j] in the innermost loop
        //       would also work

        double out = 0;
        for(k = 0; k < NB; ++k) { //NB^3 times
          out = out + A[i*NB+k] * B[k*NB+j]; //1 FMUL, 1 FADD
        }
        C[i*NB+j] += out; //1 FADD
      }
    }
} //Total opcount: 2*NB^3 + 1*NB^2

