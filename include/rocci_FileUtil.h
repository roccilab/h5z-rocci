#ifndef _ROCCI_FILE_UTIL
#define _ROCCI_FILE_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

void readfile(const char *file, const size_t num, size_t elem_size, void *data);

#endif