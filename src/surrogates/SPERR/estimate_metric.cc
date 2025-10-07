#include <SPERR/SZInterp.hpp>

double estimate_cr(double eb, float* data, double samplingRate, int32_t blockSize, std::vector<size_t> input_dims, bool skip_outliers) 
{
    if(input_dims.size() == 2) {
        return estimateSPERRCRbasedonErrorBound<float, 2>(eb, data, samplingRate, blockSize, input_dims, skip_outliers);
    }
    else { // N = 3
        return estimateSPERRCRbasedonErrorBound<float, 3>(eb, data, samplingRate, blockSize, input_dims, skip_outliers);
    }
}   

double estimate_psnr(double eb, float* data, double samplingRate, int32_t blockSize, std::vector<size_t> input_dims, bool skip_outliers) 
{
    if(input_dims.size() == 2) {
        return estimateSPERRPSNRbasedonErrorBound<float, 2>(eb, data, samplingRate, blockSize, input_dims, skip_outliers);
    }
    else { // N = 3
        return estimateSPERRPSNRbasedonErrorBound<float, 3>(eb, data, samplingRate, blockSize, input_dims, skip_outliers);
    }
}  

double estimate_ssim(double eb, float* data, double samplingRate, int32_t blockSize, std::vector<size_t> input_dims, bool skip_outliers) 
{
    if(input_dims.size() == 2) {
        return estimateSPERRSSIMbasedonErrorBound<float, 2>(eb, data, samplingRate, blockSize, input_dims, skip_outliers);
    }
    else { // N = 3
        return estimateSPERRSSIMbasedonErrorBound<float, 3>(eb, data, samplingRate, blockSize, input_dims, skip_outliers);
    }
}  