
#include <SZ3/predictor/Predictor.hpp>
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
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
class SZInterpolationCREstimator {
public:
    SZInterpolationCREstimator(Quantizer quantizer, Encoder encoder, Lossless lossless) :
            quantizer(quantizer), encoder(encoder), lossless(lossless) {

        static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quatizer interface");
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        // static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
        //               "must implement the lossless interface");
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
            sample_stride = 1;
            int n_samples = 5000;
            cmpr_err.clear();
            cmpr_err.reserve(n_samples);
            for (uint level = interpolation_level; level > 1 && level <= interpolation_level; level--) {
                size_t stride = 1U << (level - 1);
                block_interpolation(data, begin_idx, end_idx, PB_predict_collect_err,
                                    interpolators[conf.interpAlgo], conf.interpDirection, stride);
                if (cmpr_err.size() > n_samples) {
                    break;
                }
                cmpr_err.clear();
            }
        }

        {
            // the sampling process is done only on the last level, because it covers 87.5% of data points for 3D (and 75% for 2D).
            // sample_stride controls the distance of the data points covered in the sampling process.
            // original data points are used during sampling, to simulate the error impact/to make them as decompressed data,
            // errors are randomly select from cmpr_err and added for interpolation calculation.
            // TODO
            // Because sample_stride is used, the sampled quant_inds may not have same CR as the original quant_inds.
            // One possible solution is to add some zeros manually to simulate the original quant_inds
            gen = std::mt19937(rd());
            dist = std::uniform_int_distribution<>(0, cmpr_err.size());
            // sample_stride = 5;
            sample_stride = input_stride;
            // sample_stride = 2;
            quant_inds.clear();
            quant_inds.reserve(conf.num);
            block_interpolation(data, begin_idx, end_idx, PB_predict,
                                interpolators[conf.interpAlgo], conf.interpDirection, 1);
        }

        double sampling_dur = sample_timer.stop("sampling");

        encoder.preprocess_encode(quant_inds, 0);
        size_t bufferSize = 1.2 * (quantizer.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size());

        uchar *buffer = new uchar[bufferSize];
        uchar *buffer_pos = buffer;

        quantizer.save(buffer_pos);
        quantizer.postcompress_data();
        const uchar * pos = buffer_pos;

        encoder.save(buffer_pos);
        const uchar* temp = buffer_pos;
        size_t tree_size = buffer_pos - pos;
        // printf("Tree Size: %i\n", tree_size);
        // printf("QuantIndSize: %i\n", quant_inds.size());
        
        size_t huff_size = encoder.encode(quant_inds, buffer_pos);
        // printf("HuffSize: %i\n", huff_size);
        // const uchar* temp = buffer_pos-huff_size;
        auto encoder2 = SZ3::CustomHuffmanEncoder<int>();
        encoder2.load(pos, tree_size);
        auto decoded_data = encoder2.decode(temp, quant_inds.size());

        int count = 0;
        for(int i = 0; i < decoded_data.size(); i++){
            if(quant_inds[i] != decoded_data[i]){
                // printf("%i\n", i);
                count += 1;
            }
        }

        // printf("COUNT: %i\n", count);
        // printf("Unpred Huff Code: %i, %i\nZero Huff Code: %i, %i\n", encoder.get_code_for_state(0)[0], encoder.get_code_for_state(0)[1], encoder.get_code_for_state(conf.quantbinCnt / 2 )[0], encoder.get_code_for_state(conf.quantbinCnt / 2)[1]);
        // printf("Unpred Huff Code: %i\nZero Huff Code: %i\n", encoder.get_cout_for_state(0), encoder.get_cout_for_state(conf.quantbinCnt / 2));

        std::unordered_map<int,int> freq;
        int num_samples = quant_inds.size();
        int zero_bin_index = conf.quantbinCnt / 2;
        int zero_ind_freq = 0;
        int radius = 10;
        int zr_ind = (radius/2) -1; 
        int freqs[radius] = {0};
        for(int i = 0; i < num_samples; i++){
            if(quant_inds[i] == zero_bin_index){
                zero_ind_freq += 1;
                freqs[zr_ind] += 1;
            }
            else if(abs(quant_inds[i] - zero_bin_index) <= radius/2 ){
                int index;
                if(quant_inds[i] > zero_bin_index){
                    index = quant_inds[i] - zero_bin_index + zr_ind;
                }
                else if(quant_inds[i] < zero_bin_index){
                    index = zr_ind - (zero_bin_index - quant_inds[i]);
                }
                freqs[index] += 1;
            }
            freq[quant_inds[i]]++;
        }

        int center_freq = 0;
        for(int i = 0; i < radius; i++){
            int ind = i - zr_ind + zero_bin_index;
            // printf("Bin [%i]: %i | %f - %i\n", i, freqs[i], (float)freqs[i]/num_samples, encoder.get_cout_for_state(ind));
            center_freq += freqs[i];
        }

        int zero_code_len = encoder.get_cout_for_state(zero_bin_index);

        // printf("Central Bins Freq: %i | %f\n", center_freq, center_freq * 1.0 / quant_inds.size());

        std::vector<std::pair<int, int>> keyFrequencyPairs;
        int k = 5;

        // Copy key-value pairs from the unordered map to the vector
        for (const auto& kv : freq) {
            keyFrequencyPairs.push_back(std::make_pair(kv.first, kv.second));
        }
        std::sort(keyFrequencyPairs.begin(), keyFrequencyPairs.end(),
              [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
                  return p1.second > p2.second;
              });

        std::vector<int> topk_quantcodes;
        std::vector<int> topk_codelens;

        // std::cout << "Top " << k << " most frequent keys and their frequencies:" << std::endl;
        for (int i = 0; i < k && i < keyFrequencyPairs.size(); ++i) {
            // std::cout << "Key: " << keyFrequencyPairs[i].first << ", Frequency: "
            //         << keyFrequencyPairs[i].second << ", " << keyFrequencyPairs[i].second * 1.0 / num_samples
            //         << ", CodeLen: " << std::to_string(encoder.get_cout_for_state(keyFrequencyPairs[i].first)) << std::endl;

            // dont duplicate zero code estimate
            if(keyFrequencyPairs[i].first != conf.quantbinCnt /2){
                topk_quantcodes.push_back(keyFrequencyPairs[i].first);
                topk_codelens.push_back(encoder.get_cout_for_state(keyFrequencyPairs[i].first));
            }

        }

        encoder.postprocess_encode();
//            timer.stop("Coding");
        assert(buffer_pos - buffer < bufferSize);

        // printf("Estimator huffman outsize %i\n", huff_size);
		size_t postHuffmanBuffSize = buffer_pos - buffer;

        size_t compressed_size = 0;
        // uchar *lossless_data = lossless.compress(buffer,
        //                                          buffer_pos - buffer,
        //                                          compressed_size);
        uchar *lossless_data = new uchar[buffer_pos - buffer];
        compressed_size = lossless.compress(buffer, 
                                            buffer_pos - buffer, 
                                            lossless_data, 
                                            buffer_pos - buffer);
        // lossless.postcompress_data(buffer);

        // printf("Precompressedsize %i\n", compressed_size);

		// printf("Estimator huffman buffer length: %i\n", postHuffmanBuffSize);

        double prediction = 0.0;
        for (int i = 1; i < conf.quantbinCnt; i++) {
            if (freq[i] != 0) {
                float temp_bit = -log2((float)freq[i]/num_samples);
                //printf("%f %d\n", temp_bit, i);
                if (temp_bit < 32) {
                    if (temp_bit < 1) {
//                            printf("layer: %d %f\n", i, temp_bit);
                        if (i == zero_bin_index) prediction += ((float)freq[i]/num_samples) * 1;
                        else if (i == zero_bin_index-1) prediction += ((float)freq[i]/num_samples) * 2.5;
                        else if (i == zero_bin_index+1) prediction += ((float)freq[i]/num_samples) * 2.5;
                        else prediction += ((float)freq[i]/num_samples) * 4;
                    }
                    else
                        prediction += ((float)freq[i]/num_samples) * temp_bit;
                }
            }
        }
        if (freq[0] != 0) 
            prediction += ((float)freq[0]/num_samples) * 32;
            // account for uncompressed values
            compressed_size += freq[0] * sizeof(T);
            // printf("qind[0]: %i, addition: %i, %f\n", freq[0], freq[0]*sizeof(T), ((float)freq[0]/num_samples) * 32);


        float p_0 = (float) zero_ind_freq / num_samples; // percent quant inds that are zero
        float P_0; // percent of post huffman buffer size that is made up of the symbol corresponding to the zero code
        if(zero_ind_freq > num_samples/2){
            P_0 = p_0 * 1.0 / prediction;//postHuffmanBuffSize;
        }
        else {
            P_0 = -(((float)p_0) * log2((float)p_0))/ prediction;//postHuffmanBuffSize;
        }

        // printf("EST -- Zero bin: %i, Zero freq: %i, n_samples: %i, p_0: %f\n", zero_bin_index, zero_ind_freq, num_samples, p_0);

        // printf("Postcompressedsize %i\n", compressed_size);
        float pred = 1 / ((1 - p_0) * P_0 + (1 - P_0));
        float reduction_efficiency = ((1 - p_0) * P_0);
        float sum_probs = P_0;
        // float pred_ratio = postHuffmanBuffSize / pred;
        // double ratio = prediction / pred;
        // printf("ratiopred: %f\n", ratio);
        // printf("compressed_size: %i\n", compressed_size);
        // printf("PredEstBuffSize: %f, PredEstHuffmanCR: %f\n", prediction * quant_inds.size(), quant_inds.size() * sizeof(T) * 1.0 / prediction);

        // return quant_inds.size() * sizeof(T) * 1.0 / compressed_size;
        // float new_size = (prediction * quant_inds.size()) * pred_ratio;//freq[zero_bin_index];
        float huffpred = quant_inds.size() * sizeof(T) * 8.0 / (prediction * quant_inds.size());
        float lossless_pred = 32 / (prediction / pred);//quant_inds.size() * sizeof(T) * 1.0 / (new_size * 8.0);

        // compile reduction estimates for top K bins
        float pred_reduction_efficiency = reduction_efficiency;
        for(int i = 0; i < topk_quantcodes.size(); i++){
            int code = topk_quantcodes[i];
            int codelen = topk_codelens[i];
            // percent frequency
            float p_si = freq[code] * 1.0 / num_samples;
            // percent of space in huffman buffer taken by code
            float P_si = -(((float)p_si) * log2((float)p_si))/ prediction;
            // to compare empirical codelen estimation to theoretical codelen above
            float emp_P_si = ((float)p_si) * ((float)codelen) / prediction;

            pred_reduction_efficiency += ((1 - p_si) * P_si);
            sum_probs += P_si;
        }

        float pred_red_ratio = 1 / (pred_reduction_efficiency + (1-sum_probs));
        float lossless_red_ratio = 32 / (prediction / pred_red_ratio);

        // printf("newsize: %f\n", new_size);
        // printf("huffcr_pred: %f, lossless_pred: %f\n", huffpred, lossless_pred);
        // printf("lossless pred: %f, pred_red_ratio: %f\n", lossless_pred, lossless_red_ratio);

        return {huffpred, lossless_pred, sampling_dur};//lossless_pred; //32.0 / prediction;

	}
private:

    enum PredictorBehavior {
        PB_predict_collect_err, PB_predict, PB_recover
    };

    std::uniform_int_distribution<> dist;
    std::random_device rd;
    std::mt19937 gen;

    std::vector<T> cmpr_err;
    int sample_stride = 5;

    //quantize and record the quantization bins
    inline void quantize(size_t idx, T &d, T pred) {
        quant_inds.push_back(quantizer.quantize(d, pred));
    }


    //quantize and record compression error
    inline void quantize2(size_t idx, T &d, T pred) {
        T d0 = d;
        quantizer.quantize(d, pred);
        cmpr_err.push_back(d0 - d);
        d = d0;
    }

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

        return 0;
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
};

template<class T, uint N>
double estimate_compress(Config conf, T *data, double abs, int stride) {
    conf.cmprAlgo = ALGO_INTERP;
    conf.interpAlgo = INTERP_ALGO_CUBIC;
    conf.interpDirection = 0;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = abs;

	int numBins = conf.quantbinCnt / 2;

    T* data_copy = new T[conf.num];
    memcpy(data_copy, data, conf.num);

    SZInterpolationCREstimator<T, N, SZ3::LinearQuantizer<T>, SZ3::CustomHuffmanEncoder<int>, SZ3::Lossless_bypass> estimator(
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