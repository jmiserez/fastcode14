#include <stdlib.h>
#include <stdio.h>

#include "perfmon_wrapper.h"

#define COST_MEASURE // this goes before the include!
#include "cost_model.h"

COST_VARIABLES_HERE

void comp_init();
void comp_compute();
void comp_cleanup();

int main(int argc, char *argv[])
{
    int ret;

    // For a list of valid performance events, see
    //   libpfm/examples/showevtinfo -L
    // or
    //   Intel 64 and IA-32 Architectures Software Developer’s Manual,
    //     Volume 3B: System Programming Guide, Part 2
    
    // initialize perfmon
    char *events[] = {
        // on Haswell, can use at most 4, otherwise '0' results
        "PERF_COUNT_HW_CPU_CYCLES",
        "CACHE-REFERENCES",
        "CACHE-MISSES",
        "DTLB-LOAD-MISSES",
        NULL
    };

    struct perf_data *data;
    
    ret = perf_init(events, &data);
    if (ret < 0) {
        fprintf(stderr, "Could not initialize perfmon.\n");
        perf_cleanup(data);
        return 0;
    }

    // measure performance
    comp_init();
    perf_start(data);
    comp_compute();
    ret = perf_stop(data);
    if (ret < 0) {
        fprintf(stderr, "Error measuring performance.\n");
        perf_cleanup(data);
        return 0;
    }

    // report results
    printf("Performance Counters:\n");
    for (struct perf_data *iter = data; iter->name; ++iter) {
        printf("  %-50s : %"PRId64"\n", iter->name, iter->value);
    }

    printf("Cost Model:\n");
    printf("  Add: %"PRI_COST"\n", COST_ADD);
    printf("  Mul: %"PRI_COST"\n", COST_MUL);
    printf("  Div: %"PRI_COST"\n", COST_DIV);
    printf("  Pow: %"PRI_COST"\n", COST_POW);
    printf("  Abs: %"PRI_COST"\n", COST_ABS);

    printf("Derived Results\n");
    double flops      = COST_ADD + COST_MUL + COST_DIV;
    double cycles     = data[0].value; // index as in 'events'
    double cache_load = data[1].value;
    double cache_miss = data[2].value;
    printf("  Flops           : %.0lf\n", flops);
    printf("  Performance     : %.3lf\n", flops/cycles);
    printf("  Cache Miss Rate : %.3lf\n", cache_miss/cache_load);

    // cleanup
    perf_cleanup(data);
    comp_cleanup();

    return 0;
}

// Example MMM ------------------------------------------------------------

const int N = 300;
double *A, *B, *C;

void comp_init()
{
    A = malloc(N*N*sizeof(double));
    B = malloc(N*N*sizeof(double));
    C = malloc(N*N*sizeof(double));
}

void comp_compute()
{
    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            for (int k=0; k<N; ++k) {
                COST_INC_ADD(1); COST_INC_MUL(1);
                C[i*N+j] = A[i*N+k] * B[k*N+j] + C[i*N+j];
            }
        }
    }
}

void comp_cleanup()
{
    free(A);
    free(B);
    free(C);
}
