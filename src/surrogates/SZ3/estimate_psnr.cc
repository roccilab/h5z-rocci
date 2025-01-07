
#include <SZ3/predictor/Predictor.hpp>
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Extraction.hpp"
#include <SZ3/estimate_metric.h>
#include <SZ3/CustomHuffmanEncoder.hpp>
#include <memory>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <sys/stat.h>
#include <limits.h>
#include <algorithm>

using namespace SZ3;

template<class T, uint N, class Quantizer, class Encoder, class Lossless>
class SZInterpolationPSNREstimator {
public:
    SZInterpolationPSNREstimator(Quantizer quantizer, Encoder encoder, Lossless lossless) :
            quantizer(quantizer), encoder(encoder), lossless(lossless) {

        static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                      "must implement the quatizer interface");
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    std::vector<int> get_quant_inds() {
            return quant_inds;
	}

    std::vector<double> estimate(const Config &conf, T *data, int input_stride) {


        Timer sample_timer(true);
        dimension_offsets[N - 1] = 1;
        for (int i = N - 2; i >= 0; i--) {
            dimension_offsets[i] = dimension_offsets[i + 1] * conf.dims[i + 1];
        }

        dimension_sequences = std::vector<std::array<int, N>>();
        auto sequence = std::array<int, N>();
        for (int i = 0; i < N; i++) {
            sequence[i] = i;
        }
        do {
            dimension_sequences.push_back(sequence);
        } while (std::next_permutation(sequence.begin(), sequence.end()));


        std::array<size_t, N> begin_idx, end_idx;
        for (int i = 0; i < N; i++) {
            begin_idx[i] = 0;
            end_idx[i] = conf.dims[i] - 1;
        }

        int interpolation_level = -1;
        for (int i = 0; i < N; i++) {
            if (interpolation_level < ceil(log2(conf.dims[i]))) {
                interpolation_level = (int) ceil(log2(conf.dims[i]));
            }
        }

        {
            // This code block collects compression errors from a higher level into cmpr_err
            // it is observed that most levels have similar error distribution
            sample_stride = input_stride;
            size_t n_samples = 1000;
            cmpr_err.clear();
            square_err.clear();
            cmpr_err.reserve(n_samples);
            square_err.reserve(n_samples);
            for (uint level = interpolation_level; level > 1 && level <= interpolation_level; level--) {
                size_t stride = 1U << (level - 1);
                block_interpolation(data, begin_idx, end_idx, PB_predict_collect_err,
                                    interpolators[conf.interpAlgo], conf.interpDirection, stride);
                if (cmpr_err.size() > n_samples) {
                    break;
                }
                cmpr_err.clear();
                square_err.clear();
            }
        }

        double sampling_dur = sample_timer.stop("sampling");

        // {
        //     // the sampling process is done only on the last level, because it covers 87.5% of data points for 3D (and 75% for 2D).
        //     // sample_stride controls the distance of the data points covered in the sampling process.
        //     // original data points are used during sampling, to simulate the error impact/to make them as decompressed data,
        //     // errors are randomly select from cmpr_err and added for interpolation calculation.
        //     // TODO
        //     // Because sample_stride is used, the sampled quant_inds may not have same CR as the original quant_inds.
        //     // One possible solution is to add some zeros manually to simulate the original quant_inds
        //     gen = std::mt19937(rd());
        //     dist = std::uniform_int_distribution<>(0, cmpr_err.size());
        //     // sample_stride = 5;
        //     sample_stride = input_stride;
        //     // sample_stride = 2;
        //     quant_inds.clear();
        //     quant_inds.reserve(conf.num);
        //     block_interpolation(data, begin_idx, end_idx, PB_predict,
        //                         interpolators[conf.interpAlgo], conf.interpDirection, 1);
        // }

        maxVal = data[0];
        minVal = data[0];
        for(int k = 0; k < conf.num; k++){
            if(data[k] > maxVal) maxVal = data[k];
            if(data[k] < minVal) minVal = data[k];
        }

        auto const count = static_cast<float>(square_err.size());
        double eps = 1e-16;
        double mse = std::reduce(square_err.begin(), square_err.end()) / count;
        double value_range = maxVal - minVal;
        double estPSNR = -20.0*log10((sqrt(mse)/value_range) + eps);

        double mean_err = std::reduce(cmpr_err.begin(), cmpr_err.end()) / count;

        return {estPSNR, sampling_dur};

	}
private:

    enum PredictorBehavior {
        PB_predict_collect_err, PB_predict, PB_recover
    };

    std::uniform_int_distribution<> dist;
    std::random_device rd;
    std::mt19937 gen;

    std::vector<T> cmpr_err;
    std::vector<T> square_err;
    int sample_stride = 5;

    //quantize and record the quantization bins
    inline void quantize(size_t idx, T &d, T pred) {
        quant_inds.push_back(quantizer.quantize(d, pred));
    }

    //quantize and record compression error
    inline void quantize2(size_t idx, T &d, T pred) {
        T d0 = d;
        if(d0 > maxVal) maxVal = d0;
        if(d0 < minVal) minVal = d0;
        quantizer.quantize_and_overwrite(d, pred);
        cmpr_err.push_back(d0 - d);
        square_err.push_back(pow(d0 - d, 2));
        d = d0;
    }

    //recover and return the squared error in prediction
    // inline T recover(size_t idx, T &d, T pred) {
    //     T pred = quantizer.recover(pred, quant_inds[quant_index++]);
    //     return pow(pred - d, 2);
    // };

    //Add noise/compression error to original data, to simulate it as decompressed data
    inline T s(T d) {
        return d + cmpr_err[dist(gen)];
    }

    double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,
                                  const std::string &interp_func,
                                  const PredictorBehavior pb) {
        size_t n = (end - begin) / stride + 1;
        if (n <= 1) {
            return 0;
        }

        size_t stride3x = 3 * stride;
        size_t stride5x = 5 * stride;

        double predict_error = 0;

        if (pb == PB_predict_collect_err) {
            if (interp_func == "linear" || n < 5) {
                for (size_t i = 1; i + 1 < n; i += 2 * sample_stride) {
                    T *d = data + begin + i * stride;
                    quantize2(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 4) {
                        quantize2(d - data, *d, *(d - stride));
                    } else {
                        quantize2(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                    }
                }
            } else {
                T *d;
                size_t i;
                for (i = 3; i + 3 < n; i += 2 * sample_stride) {
                    d = data + begin + i * stride;
                    quantize2(d - data, *d,
                                               interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                }
                d = data + begin + stride;
                quantize2(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));


                d = data + begin + ((n % 2 == 0) ? n - 3 : n - 2) * stride;
                quantize2(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                if (n % 2 == 0) {
                    d = data + begin + (n - 1) * stride;
                    quantize2(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                }

            }
        } else {
            if (interp_func == "linear" || n < 5) {
                for (size_t i = 1; i + 1 < n; i += 2 * sample_stride) {
                    T *d = data + begin + i * stride;
                    quantize(d - data, *d, interp_linear(s(*(d - stride)), s(*(d + stride))));
                }
                if (n % 2 == 0) {
                    T *d = data + begin + (n - 1) * stride;
                    if (n < 4) {
                        quantize(d - data, *d, s(*(d - stride)));
                    } else {
                        quantize(d - data, *d, interp_linear1(s(*(d - stride3x)), s(*(d - stride))));
                    }
                }
            } else {
                T *d;
                size_t i;
                for (i = 3; i + 3 < n; i += 2 * sample_stride) {
                    d = data + begin + i * stride;
                    quantize(d - data, *d,
                                              interp_cubic(s(*(d - stride3x)), s(*(d - stride)), s(*(d + stride)), s(*(d + stride3x))));
                }
                d = data + begin + stride;
                quantize(d - data, *d, interp_quad_1(s(*(d - stride)), s(*(d + stride)), s(*(d + stride3x))));


                d = data + begin + ((n % 2 == 0) ? n - 3 : n - 2) * stride;
                quantize(d - data, *d, interp_quad_2(s(*(d - stride3x)), s(*(d - stride)), s(*(d + stride))));
                if (n % 2 == 0) {
                    d = data + begin + (n - 1) * stride;
                    quantize(d - data, *d, interp_quad_3(s(*(d - stride5x)), s(*(d - stride3x)), s(*(d - stride))));
                }
            }
        } 

        return predict_error;
    }

    template<uint NN = N>
    typename std::enable_if<NN == 1, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
    }

    template<uint NN = N>
    typename std::enable_if<NN == 2, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x * sample_stride) {
            size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
            predict_error += block_interpolation_1d(data, begin_offset,
                                                    begin_offset +
                                                    (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                    stride * dimension_offsets[dims[0]], interp_func, pb);
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride * sample_stride) {
            size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
            predict_error += block_interpolation_1d(data, begin_offset,
                                                    begin_offset +
                                                    (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                    stride * dimension_offsets[dims[1]], interp_func, pb);
        }
        return predict_error;
    }

    template<uint NN = N>
    typename std::enable_if<NN == 3, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x * sample_stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x * sample_stride) {
                size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                      k * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[0]] - begin[dims[0]]) *
                                                        dimension_offsets[dims[0]],
                                                        stride * dimension_offsets[dims[0]], interp_func, pb);
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride * sample_stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x * sample_stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                      k * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[1]] - begin[dims[1]]) *
                                                        dimension_offsets[dims[1]],
                                                        stride * dimension_offsets[dims[1]], interp_func, pb);
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride * sample_stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride * sample_stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                      begin[dims[2]] * dimension_offsets[dims[2]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[2]] - begin[dims[2]]) *
                                                        dimension_offsets[dims[2]],
                                                        stride * dimension_offsets[dims[2]], interp_func, pb);
            }
        }
        return predict_error;
    }


    template<uint NN = N>
    typename std::enable_if<NN == 4, double>::type
    block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                        const std::string &interp_func, const int direction, size_t stride = 1) {
        double predict_error = 0;
        size_t stride2x = stride * 2;
        const std::array<int, N> dims = dimension_sequences[direction];
        for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x * sample_stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x * sample_stride) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                     t <= end[dims[3]]; t += stride2x * sample_stride) {
                    size_t begin_offset =
                            begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                            k * dimension_offsets[dims[2]] +
                            t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[0]] - begin[dims[0]]) *
                                                            dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[0]], interp_func, pb);
                }
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride * sample_stride) {
            for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x * sample_stride) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                     t <= end[dims[3]]; t += stride2x * sample_stride) {
                    size_t begin_offset =
                            i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                            k * dimension_offsets[dims[2]] +
                            t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[1]] - begin[dims[1]]) *
                                                            dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb);
                }
            }
        }
        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride * sample_stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride * sample_stride) {
                for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                     t <= end[dims[3]]; t += stride2x * sample_stride) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          begin[dims[2]] * dimension_offsets[dims[2]] +
                                          t * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[2]] - begin[dims[2]]) *
                                                            dimension_offsets[dims[2]],
                                                            stride * dimension_offsets[dims[2]], interp_func, pb);
                }
            }
        }

        for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride * sample_stride) {
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride * sample_stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride * sample_stride) {
                    size_t begin_offset =
                            i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                            k * dimension_offsets[dims[2]] +
                            begin[dims[3]] * dimension_offsets[dims[3]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[3]] - begin[dims[3]]) *
                                                            dimension_offsets[dims[3]],
                                                            stride * dimension_offsets[dims[3]], interp_func, pb);
                }
            }
        }
        return predict_error;
    }

public:
    std::vector<std::string> interpolators = {"linear", "cubic"};
    std::array<size_t, N> dimension_offsets;
    std::vector<std::array<int, N>> dimension_sequences;
    Quantizer quantizer;
    Encoder encoder;
    Lossless lossless;
    std::vector<int> quant_inds;
    float maxVal = std::numeric_limits<float>::min(); // store max value observed in data;
    float minVal = std::numeric_limits<float>::max();
};

template<class T, uint N>
double estimate_psnr(Config conf, T *data, double abs, int stride) {
    conf.cmprAlgo = ALGO_INTERP;
    conf.interpAlgo = INTERP_ALGO_CUBIC;
    conf.interpDirection = 0;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs;

	int numBins = conf.quantbinCnt / 2;

    T* data_copy = new T[conf.num];
    memcpy(data_copy, data, sizeof(T)*conf.num);

    SZInterpolationPSNREstimator<T, N, SZ3::LinearQuantizer<T>, SZ3::CustomHuffmanEncoder<int>, SZ3::Lossless_bypass> estimator(
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