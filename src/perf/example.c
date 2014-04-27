#include <stdlib.h>
#include <stdio.h>

#include "perfmon_wrapper.h"

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
	for (struct perf_data *iter = data; iter->name; ++iter) {
		printf("%-50s : %"PRId64"\n", iter->name, iter->value);
	}

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
