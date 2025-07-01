#include <cstdlib>
#include <cstring>

struct FloatPtrArray {
    float *arr[2];
};

double SZp_float_estimate_compress(float* data, size_t nbEle, int blockSize, float absErrBound, size_t sample_stride);
double SZp_float_estimate_psnr(float* data, size_t nbEle, int blockSize, float absErrBound, size_t sample_stride);
struct FloatPtrArray SZp_float_estimate_ssim(float* data, size_t dims[3], size_t nbEle, size_t n_dim, int blockSize, float absErrBound, size_t sample_stride);
