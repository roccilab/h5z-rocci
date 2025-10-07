#include <sys/time.h>

float szx_estimate_cr_float(const float* data, float absErrorBound, int samplingStride, size_t nbEle);
float szx_estimate_ssim_float(const float* data, float absErrorBound, int samplingStride, int blockSize, size_t dims[3], size_t nbEle);
float szx_estimate_psnr_float(const float* data, float absErrorBound, int samplingStride, size_t nbEle);
