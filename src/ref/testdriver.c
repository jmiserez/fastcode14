#include<stdio.h>
#include<stdlib.h>

typedef struct {
  char* (*f)();
  char* desc;
} testfunc_t;

// include list of files

// #include "matrix_test.h"

/*
 * Returning an empty string means success. Any non-empty string means failure.
 * In case of a failure, the test function is supposed to return a string
 * describing the reason of failure.
 *
 */
char* dummy_success() {
    return "";
}

char* dummy_fail() {
    return "Some error description;";
}

char* dummy_null() {
    return NULL;
}

size_t register_test( testfunc_t** tests, size_t n, char* (*f)(), char* desc ) {
    size_t new_n = n+1;
    *tests = (testfunc_t*) realloc( *tests, (new_n) * sizeof( testfunc_t ) );
    (*tests)[new_n-1] = (testfunc_t) { f, desc };
    return new_n;
}

// register the function that you want to test here!

size_t register_tests( testfunc_t** tests, size_t n ) {
    n = register_test( tests, n, dummy_success, "dummy_success" );
    n = register_test( tests, n, dummy_fail, "dummy_fail" );
    n = register_test( tests, n, dummy_null, "dummy_null" );
}

void run_tests( testfunc_t* tests, size_t n ) {
    size_t i;
    for( i = 0; i < n; i++ ) {
        char* res = tests[i].f();
        if( res != NULL )
            if( *res != '\0' ) {
                printf("%s: FAILED\n", tests[i].desc);
                printf("  Reason: %s\n", res);
            } else
                printf("%s: SUCCEEDED\n", tests[i].desc);
        else
            printf("warning: %s returned null!\n", tests[i].desc);
    }
}

int main( int argc, char* argv[] ) {
    testfunc_t* tests = NULL;
    size_t test_count = 0;

    test_count = register_tests( &tests, test_count );
    run_tests( tests, test_count );

    free(tests);
    return 0;
}
