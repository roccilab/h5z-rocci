#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

void sample_1d_float(float* data, size_t blocksize, size_t sample_stride, size_t dims[], float** outBuffer, size_t* outSize, size_t sample_dims[]);
void sample_2d_float(float* data, size_t blocksize, size_t sample_stride, size_t dims[], float** outBuffer, size_t* outSize, size_t sample_dims[]);
void sample_3d_float(float* data, size_t blocksize, size_t sample_stride, size_t dims[], float** outBuffer, size_t* outSize, size_t sample_dims[]);