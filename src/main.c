#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "rdtsc.h"

//defined from assignment
#ifndef NB
#error NB is not defined
#endif

//defined in Makefile
#ifndef CPUMHZ
#error CPUMHZ is not defined
#endif

//definitions from hw02, make sure each measurement is at least 1 second
#define RUNS     2
#define CYCLES_REQUIRED (CPUMHZ*1000000)


double *A, *B, *C, *Cref;

void compute();

void fill(){
    int i;
#ifndef N
    for(i = 0; i<NB*NB; i++) {
#else
    for(i = 0; i<N*N; i++) {
#endif
      A[i] = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
      B[i] = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
      C[i] = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
      Cref[i] = C[i];
    }
}

void verify() {
  int i, j, k;
#ifndef N
  for(i = 0; i < NB; ++i) {
    for(j = 0; j < NB; ++j) {
      double out = 0;
      for(k = 0; k < NB; ++k) {
        out = out + A[i*NB+k] * B[k*NB+j];
      }
      Cref[i*NB+j] += out;
    }
  }
  for(i = 0; i<NB*NB; i++) {
    if (abs(Cref[i]-C[i])>10e-7) {
      printf("error\n");
      exit(EXIT_FAILURE);
    }
  }
#else
  for(i = 0; i < N; ++i) {
    for(j = 0; j < N; ++j) {
      double out = 0;
      for(k = 0; k < N; ++k) {
        out = out + A[i*N+k] * B[k*N+j];
      }
      Cref[i*N+j] += out;
    }
  }
  for(i = 0; i<N*N; i++) {
    if (abs(Cref[i]-C[i])>10e-7) {
      printf("error\n");
      exit(EXIT_FAILURE);
    }
  }
#endif
}

/**
 * Uses the RDTSC method from hw02. This method was very precise well on my machine in hw02.
 */
void benchmark(){
    tsc_counter start, end;
    double cycles = 0.;
    size_t num_runs = RUNS, multiplier = 1;
    int i;


    //Cache warm-up + timing setup
    do{
        fill();
        compute();
        num_runs = num_runs * multiplier;
        CPUID(); RDTSC(start);
        for(i = 0; i < num_runs; i++) {
          compute();
        }
        CPUID(); RDTSC(end);

        cycles = (double)(COUNTER_DIFF(end, start));
        multiplier = ceil (  (CYCLES_REQUIRED) / (cycles)  + 1.0 );

    }while (multiplier > 2);

    fill();
    compute();
    compute();
    compute();
    compute();
    compute();
    CPUID(); RDTSC(start);
    for (i = 0; i < num_runs; ++i) {
      compute();
    }
    CPUID(); RDTSC(end);

    cycles = (double)(COUNTER_DIFF(end, start))/num_runs;

    printf("\nnum_runs:       %lu\n", num_runs);
#ifdef N
    printf("N:                %d\n", N);
#endif
    printf("NB:               %d\n", NB);
    printf("Speed [MHz]:      %.2f\n", CPUMHZ);
    printf("Op count [flops]: %d\n", OPCOUNT);
    printf("Runtime [cycles]: %.2f\n", cycles);
    printf("Benchmark time [s]: %.2f\n", (double)(COUNTER_DIFF(end, start))/(((double)CPUMHZ)*(double)1000000.0));
    printf("Performance [flops/cycle]: %.2f\n", ((double)OPCOUNT)/cycles);
}

int main(){
  int i;
#ifndef N
  A = (double*)malloc(NB*NB*sizeof(double));
  B = (double*)malloc(NB*NB*sizeof(double));
  C = (double*)malloc(NB*NB*sizeof(double));
  Cref = (double*)malloc(NB*NB*sizeof(double));
#else
  A = (double*)malloc(N*N*sizeof(double));
  B = (double*)malloc(N*N*sizeof(double));
  C = (double*)malloc(N*N*sizeof(double));
  Cref = (double*)malloc(N*N*sizeof(double));
#endif

  srand ( time(NULL) );

  //verify 5 random matrices
  for(i = 0; i < 5; i++){
      fill();
      compute();
      verify();
  }
  printf("Verification successful!\n");
  benchmark();

  return 0;
}
