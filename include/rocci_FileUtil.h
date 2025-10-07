#ifndef _ROCCI_FILE_UTIL
#define _ROCCI_FILE_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

void readfile(const char *file, const size_t num, size_t elem_size, void *data);

#ifdef __cplusplus
}
#endif
#endif