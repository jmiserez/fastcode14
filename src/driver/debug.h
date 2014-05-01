#include <stdio.h>

#ifndef DRV_DEBUG_H
#define DRV_DEBUG_H

#define FUSION_ERR(...)                                              \
        fprintf(stderr, "%s:%d:%s: ", __FILE__, __LINE__, __func__); \
        fprintf(stderr, __VA_ARGS__);

#endif // DRV_DEBUG_H
