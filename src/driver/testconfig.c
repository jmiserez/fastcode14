#include"testconfig.h"
#include"image_io.h"
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#include<ctype.h>
#include<stdio.h>
#include<assert.h>

/* Yay! String parsing in C! <3 */

char* strncpy_nt( char* dest, const char* p_start, const char* p_stop ) {
    assert(p_stop >= p_start);
    size_t len = p_stop - p_start;
    if( dest ) {
        strncpy(dest, p_start, len);
        dest[len] = '\0';
    }
    return dest;
}

char* strncpy_alloc( const char* p_start, const char* p_stop ) {
    assert(p_stop >= p_start);
    size_t len = p_stop - p_start;
    char* dst;

    dst = (char*) malloc( len * sizeof(char) + 1 );
    if( dst )
        strncpy_nt(dst, p_start, p_stop);
    return dst;
}

int readline( char* buf, size_t n, FILE* f ) {
    int length = 0;
    char ch = getc(f);
    while ( (ch != '\n') && (ch != EOF) && (length < (n-1)) ) {
        if( !isspace(ch) ) { // ignore whitespace
            buf[length] = ch;
            length++;
        }
        ch = getc(f);
    }
    buf[length] = '\0';
    return length;
}

char* read_fparam( double* param, char* p_start ) {
    const int max_double_len = 18;
    char double_str[max_double_len+1];
    char* p_cur = p_start;
    int point = 0;

    *param = 0.0;

    while( (*p_cur != '\0')
           && (isdigit(*p_cur) || (*p_cur=='.' && !point))
           && (p_cur-p_start < max_double_len)) {
        if( *p_cur == '.' ) point = 1;
        p_cur++;
    }

    strncpy_nt( double_str, p_start, p_cur );
    *param = atof( double_str );

    return p_cur;
}

int parse_filename( testconfig_t* testconfig, char* filename ) {
    int i;
    char* p_cur = filename;
    char* p_start = filename;

    const int max_int_len = 9;
    char n_files_str[max_int_len+1];

    double *params[3] = {&testconfig->contrast, &testconfig->saturation, &testconfig->exposure};

    /* read the prefix */
    while( (*p_cur != '\0') && (*p_cur != '-') )
        p_cur++;

    char* prefix = strncpy_alloc(p_start, p_cur);
    testconfig->prefix = prefix;
    // testconfig->prefix

    // parse error if prefix is not followed by '-'
    if( *p_cur != '-' ) return (int) (filename - p_cur);
    p_cur++;

    /* # read number_of_files # */

    p_start = p_cur;
    while( (*p_cur != '\0') && (*p_cur != '-') && (isdigit(*p_cur))
           && ((p_cur - p_start) < max_int_len))
        p_cur++;
    // found non-digit character => parse_error
    if( (*p_cur != '-') && !isdigit(*p_cur) ) return (int) (filename - p_cur);
    strncpy_nt( n_files_str, p_start, p_cur );
    testconfig->number_of_files = atoi( n_files_str );
    p_cur++;

    /* # read float parameters */

    // contrast and saturation
    for( i = 0; i < 2; i++ ) {
        p_start = p_cur;
        p_cur = read_fparam( params[i], p_start );
        if( (*p_cur != '-') || (p_cur == p_start) ) return (int) (filename - p_cur);
        p_cur++;
    }

    // exposure
    p_start = p_cur;
    p_cur = read_fparam( &testconfig->exposure, p_start );
    if( p_cur == p_start ) return (int) (filename - p_cur);

    // if p_cur points to a dot, we stand on the extension separator
    if( *p_cur == '.' ) p_cur++;

    p_start = p_cur;
    while( *p_cur != '\0' )
        p_cur++;

    testconfig->extension = strncpy_alloc( p_start, p_cur );
    testconfig->filename = strncpy_alloc( filename, p_cur );

    return p_cur-filename;
}

int read_testconfigurations( testconfig_t* testconfigs, size_t max_configs, FILE* f, char* srcPath, char* refPath ) {
    int config_count = 0;
    int read_line_sz = 0;
    int ret_parse_fn;
    testconfig_t* testconfig = testconfigs;
    const int line_sz = 256;
    char line[line_sz];
    char line_count = 0;
    size_t max_len = 0;
    int i;

    if( f != NULL ) {
        while( !feof(f) && (config_count < max_configs) ) {
            line_count++;
            if( ((read_line_sz = readline(line, line_sz, f)) > 0)
                    && (line[0] != '#') ) { // exclude comments
                if( (ret_parse_fn = parse_filename( testconfig, line )) > 0 ) {

                    // for convenience, we set the full paths here
                    max_len = strlen(refPath) + strlen(testconfig->filename) + 2; // plus '/' and '/0'
                    testconfig->ref_path = (char*) malloc( max_len * sizeof(char) );
                    snprintf(testconfig->ref_path, max_len, "%s/%s", refPath, testconfig->filename);

                    // note that the input file names cannot be longer than the reference filenames.
                    max_len = strlen(srcPath) + strlen(testconfig->filename) + 2; // plus '/' and '/0'
                    testconfig->input_paths = (char**) malloc( testconfig->number_of_files * sizeof(char*) );
                    for( i = 0; i < testconfig->number_of_files; i++ ) {
                        testconfig->input_paths[i] = (char*) malloc( max_len * sizeof( char* ) );
                        snprintf(testconfig->input_paths[i], max_len, "%s/%s.%d.%s", srcPath,
                                 testconfig->prefix, i, testconfig->extension);
                    }

                    testconfig++;
                    config_count++;
                } else {
#ifdef DEBUG
                    fprintf(stderr, "error at line %d col %d parsing line:\n%s\n",
                            line_count, -ret_parse_fn, line);
#endif
                }
            }
        }
    }
    return config_count;
}

void print_testconfiguration( testconfig_t* tc ) {
    int i;
    printf("filename:     %s\n", tc->filename);
    printf("prefix:       %s\n", tc->prefix);
    printf("extension:    %s\n", tc->extension);
    printf("no. of files: %ld\n", tc->number_of_files);
    printf("contrast:     %lf\n", tc->contrast);
    printf("saturation:   %lf\n", tc->saturation);
    printf("exposure:     %lf\n", tc->exposure);
    printf("ref. path:    %s\n", tc->ref_path);
    printf("input paths:\n");
    for( i = 0; i < tc->number_of_files; i++ ) {
        printf("              %s\n", tc->input_paths[i]);
    }
}

double** tc_read_input_images( size_t* read_imgs, uint32_t *ret_w, uint32_t *ret_h,
                               testconfig_t* tc ) {
    int i;
    uint32_t w, h, new_w, new_h;

    double** images = (double**) malloc( tc->number_of_files * sizeof(double*));

    if( tc->number_of_files > 0 ) {
        images[0] = load_tiff_rgb( &w, &h, tc->input_paths[0] );
    }
    if( images[0] != NULL ) {
        for( i = 1; i < tc->number_of_files; i++ ) {
            images[i] = load_tiff_rgb( &new_w, &new_h, tc->input_paths[i] );
            if( !(new_w == w && new_h == h) || (images[i] == NULL) ) {
                free(images[i]);
                i--;
                break;
            }
            w = new_w; h = new_h;
        }
    } else
        free( images[0] );

    *ret_w = w; *ret_h = h;

    if( i < tc->number_of_files )
        images = realloc( images, i * sizeof(double*));
    *read_imgs = i;

    return images;
}

void tc_free_input_images( double** images, size_t n_images ) {
    int i;
    for( i = 0; i < n_images; i++ ) {
        free_rgb( images[i] );
    }
    free(images);
}

void tc_free( testconfig_t* tc ) {
    int i;
    free(tc->filename);
    free(tc->prefix);
    free(tc->extension);
    free(tc->ref_path);
    for( i = 0; i < tc->number_of_files; i++ ) {
        free(tc->input_paths[i]);
    }
    free(tc->input_paths);
}

void tcs_free( testconfig_t* tc, size_t n_tc ) {
    size_t i;
    for( i = 0; i < n_tc; i++ ) {
        tc_free(&tc[i]);
    }
}
