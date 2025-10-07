#include <cstdlib>
#include <vector>

double estimate_cr(double eb, float* data, double samplingRate, int32_t blockSize, std::vector<size_t> input_dims, bool skip_outliers);
double estimate_psnr(double eb, float* data, double samplingRate, int32_t blockSize, std::vector<size_t> input_dims, bool skip_outliers);
double estimate_ssim(double eb, float* data, double samplingRate, int32_t blockSize, std::vector<size_t> input_dims, bool skip_outliers);
