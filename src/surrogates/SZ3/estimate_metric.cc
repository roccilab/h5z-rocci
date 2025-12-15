#include <SZ3/estimate_metric.h>

#include "estimate_cr.cc"
#include "estimate_psnr.cc"
#include "estimate_ssim.cc"

using namespace SZ3;

template<class T, uint N>
double estimate_compress(Config conf, T *data, double abs, int stride) {
    conf.cmprAlgo = ALGO_INTERP;
    conf.interpAlgo = INTERP_ALGO_CUBIC;
    conf.interpDirection = 0;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs;

	int numBins = conf.quantbinCnt / 2;

    T* data_copy = new T[conf.num];
    memcpy(data_copy, data, conf.num * sizeof(T));

    SZ3_CR::SZInterpolationCREstimator<T, N, SZ3::LinearQuantizer<T>, SZ3::CustomHuffmanEncoder<int>, SZ3::Lossless_bypass> estimator(
            SZ3::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
            SZ3::CustomHuffmanEncoder<int>(),
            SZ3::Lossless_bypass());
    std::vector<double> estresults = estimator.estimate(conf, data, stride);
    double estCR = estresults[0];
    delete[] data_copy;

    return estCR;
}

double estimate_cr_float(SZ3::Config conf, float *data, double abs, int stride, uint N) {
    switch(N) {
        case 1:
            return estimate_compress<float, 1>(conf, data, abs, stride);
        case 2:
            return estimate_compress<float, 2>(conf, data, abs, stride);
        case 3:
            return estimate_compress<float, 3>(conf, data, abs, stride);
        case 4:
            return estimate_compress<float, 4>(conf, data, abs, stride);
        default:
            throw std::invalid_argument("N must be in range [1,4]");
    }
}


template<class T, uint N>
double estimate_psnr(Config conf, T *data, double abs, int stride) {

    printf("stride = %d\n", stride);
    fflush(stdout);
    conf.cmprAlgo = ALGO_INTERP;
    conf.interpAlgo = INTERP_ALGO_CUBIC;
    conf.interpDirection = 0;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs;

	int numBins = conf.quantbinCnt / 2;

    T* data_copy = new T[conf.num];
    memcpy(data_copy, data, sizeof(T)*conf.num);

    SZ3_PSNR::SZInterpolationPSNREstimator<T, N, SZ3::LinearQuantizer<T>, SZ3::CustomHuffmanEncoder<int>, SZ3::Lossless_bypass> estimator(
            SZ3::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
            SZ3::CustomHuffmanEncoder<int>(),
            SZ3::Lossless_bypass());
    std::vector<double> estresults = estimator.estimate(conf, data_copy, stride);
    double estPSNR = estresults[0];
    double sample_dur = estresults[1];
    delete[] data_copy;
	
    return estPSNR;
}

double estimate_psnr_float(SZ3::Config conf, float *data, double abs, int stride, uint N) {
    switch(N) {
        case 1:
            return estimate_psnr<float, 1>(conf, data, abs, stride);
        case 2:
            return estimate_psnr<float, 2>(conf, data, abs, stride);
        case 3:
            return estimate_psnr<float, 3>(conf, data, abs, stride);
        case 4:
            return estimate_psnr<float, 4>(conf, data, abs, stride);
        default:
            throw std::invalid_argument("N must be in range [1,4]");
    }
}

template<class T, uint N>
double estimate_ssim(Config conf, T *data, double abs, int stride, int blocksize) {
    conf.cmprAlgo = ALGO_INTERP;
    conf.interpAlgo = INTERP_ALGO_CUBIC;
    conf.interpDirection = 0;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs;

	int numBins = conf.quantbinCnt / 2;

    T* data_copy = new T[conf.num];
    memcpy(data_copy, data, sizeof(T)*conf.num);

    SZ3_SSIM::SZInterpolationSSIMEstimator<T, N, SZ3::LinearQuantizer<T>, SZ3::CustomHuffmanEncoder<int>, SZ3::Lossless_bypass> estimator(
            SZ3::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
            SZ3::CustomHuffmanEncoder<int>(),
            SZ3::Lossless_bypass());
    std::vector<double> estresults = estimator.estimate(conf, data_copy, stride, blocksize);
    double estSSIM = estresults[0];
    double sample_dur = estresults[1];
    delete[] data_copy;

    return estSSIM;

}

double estimate_ssim_float(SZ3::Config conf, float *data, double abs, int stride, int blocksize, uint N) {
    switch(N) {
        case 1:
            return estimate_ssim<float, 1>(conf, data, abs, stride, blocksize);
        case 2:
            return estimate_ssim<float, 2>(conf, data, abs, stride, blocksize);
        case 3:
            return estimate_ssim<float, 3>(conf, data, abs, stride, blocksize);
        default:
            throw std::invalid_argument("N must be in range [1,3]");
    }
}