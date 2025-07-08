#ifndef __SELECT_HPP__
#define __SELECT_HPP__

#include <unordered_map>

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <type_traits>

#include <libpressio.h>
#include <libpressio_ext/cpp/libpressio.h>
#include <libpressio_ext/cpp/compressor.h>
#include <libpressio_ext/cpp/options.h>
#include <libpressio_ext/cpp/data.h>
#include <libpressio_meta.h>

#include "rocci_config.hpp"

#include "surrogates/SZ3/SZ3Surrogate.cc"
#include "surrogates/SZx/SZxSurrogate.cc"
#include "surrogates/ZFP/ZFPSurrogate.cc"
#include "surrogates/SPERR/SPERRSurrogate.cc"
#include "surrogates/SZp/SZpSurrogate.cc"
// #include "surrogates/cuSZp/cuSZpSurrogate.cc"

pressio library;

using namespace std::string_literals;

// static const auto compressor_list = {"sz3", "zfp", "szx", "sperr"};
// static const auto compressor_list = {"sz3", "zfp", "szx"};
static const auto compressor_list = {"sz3"};

std::map<std::string, uint8_t> is_reversed_surrogate_dimension = {
    {"sz3", 1}, 
    {"zfp", 0}, 
    {"szx", 0}, 
    {"sperr", 0}, 
    {"szp", 1}
};

std::map<std::string, ALGO> compressor_string_to_enum = {
    {"sz3", ALGO::ALGO_SZ3}, 
    {"zfp", ALGO::ALGO_ZFP}, 
    {"szx", ALGO::ALGO_SZx}, 
    {"sperr", ALGO::ALGO_SPERR}, 
    {"szp", ALGO::ALGO_SZp}
};

std::map<std::string, uint8_t> &is_reversed_dimension = is_reversed_surrogate_dimension;

std::map<std::string, uint8_t> is_reversed_cmp_dimension = {
    {"sz3", 0}, 
    {"zfp", 0}, 
    {"szx", 0}, 
    {"sperr", 0}, 
    {"szp", 0}
};

double get_surrogate_cr(const std::string _, pressio_data *input, const double eb);
double get_real_cr(const std::string _, pressio_data *input, const double eb);
size_t compress_to_destination(const std::string _, pressio_data *input, const double eb, char *&output);

size_t compress_fix_eb_best_cr(pressio_data *input, const double eb, char *&output);
// size_t compress_fix_psnr_best_cr(pressio_data *input, const double psnr, char *output);
// size_t compress_fix_ssim_best_cr(pressio_data *input, const double ssim, char *output);

// size_t compress_fix_cr_best_eb(pressio_data *input, const double cr, char *output);
// size_t compress_fix_cr_best_psnr(pressio_data *input, const double cr, char *output);
// size_t compress_fix_cr_best_ssim(pressio_data *input, const double cr, char *output);

/*
./tools/select \
/anvil/projects/x-cis240669/kai/data/hurricane-100x500x500/Uf48.bin.dat \
3 500 500 100 \
~/scratch/tem eb_cr \
1e-1
*/

using function_t = std::function<size_t(
    pressio_data *, 
    const double, 
    char *&)>;

std::unordered_map<std::string, function_t> func_mp = {
    {"eb_cr", compress_fix_eb_best_cr}, 
    // {"psnr_cr", compress_fix_psnr_best_cr}, 
    // {"ssim_cr", compress_fix_ssim_best_cr}, 
    // {"cr_eb", compress_fix_cr_best_eb}, 
    // {"cr_psnr", compress_fix_cr_best_psnr}, 
    // {"cr_ssim", compress_fix_cr_best_ssim}
};

std::vector<size_t> dims, reversed_dims;

template <typename T>
size_t ROCCI_select_compress(Config &conf, T *data, char *&compressedData, size_t &compressed_size) {

    dims = conf.dims;
    reversed_dims = dims;
    std::reverse(reversed_dims.begin(), reversed_dims.end());

    double compression_parameter = conf.absErrorBound;

    pressio_dtype type = pressio_float_dtype;
    {
        if constexpr(std::is_same_v<T, float>) type = pressio_float_dtype;
        else if constexpr(std::is_same_v<T, double>) type = pressio_double_dtype;
        else if constexpr(std::is_same_v<T, std::int8_t>) type = pressio_int8_dtype;
        else if constexpr(std::is_same_v<T, std::int16_t>) type = pressio_int16_dtype;
        else if constexpr(std::is_same_v<T, std::int32_t>) type = pressio_int32_dtype;
        else if constexpr(std::is_same_v<T, std::int64_t>) type = pressio_int64_dtype;
        else if constexpr(std::is_same_v<T, std::uint8_t>) type = pressio_uint8_dtype;
        else if constexpr(std::is_same_v<T, std::uint16_t>) type = pressio_uint16_dtype;
        else if constexpr(std::is_same_v<T, std::uint32_t>) type = pressio_uint32_dtype;
        else if constexpr(std::is_same_v<T, std::uint64_t>) type = pressio_uint64_dtype;
        else {
            throw std::runtime_error("Unsupported data type");
        }
    }

    pressio_data input = pressio_data::nonowning(type, data, dims);
    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

    return compress_fix_eb_best_cr(&input, compression_parameter, compressedData);
}

size_t compress_fix_eb_best_cr(pressio_data *input, const double eb, char *&output) {

    double best_cr = 0;
    std::string best_compressor = "null";
    // test the best compressor
    for (const auto &_ : compressor_list) {
        auto current_cr = get_surrogate_cr(_, input, eb);
        if (current_cr > best_cr) {
            best_cr = current_cr;
            best_compressor = _;
        }
    }

    // std::cout << "Best compressor: " << best_compressor << " with CR: " << best_cr << std::endl;

    size_t compressed_size = compress_to_destination(best_compressor, input, eb, output);
    return compressed_size;
}

// auto get_surrogate(const std::string &_) {

// }

double get_surrogate_cr(const std::string _, pressio_data *input, const double eb) {

    if (is_reversed_dimension[_]) {
        pressio_data_reshape(input, reversed_dims.size(), reversed_dims.data());
    }
    else {
        pressio_data_reshape(input, dims.size(), dims.data());
    }

    const std::string &compressor_string = _;
    auto compressor = library.get_compressor(_ + "_surrogate");
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " surrogate" << std::endl;
        exit(-1);
    }

    pressio_options compression_options{
        {_ + "_surrogate:abs_error_bound", eb}, 
        // {_ + "_surrogate:error_bound_mode", "ABS"}, 
        {_ + "_surrogate:metric", "cr"}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

    if(compressor->compress(input, &compressed) != 0) {
        std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    auto metrics = compressor->get_metrics_results();
    auto assert_defined = [&](const char* key, auto value){
        if(metrics.get(key, value)!= pressio_options_key_set) {
            throw std::runtime_error(key + " was note set"s);
        }
    };

    double cr;
    assert_defined((_ + "_surrogate:cr").c_str(), &cr);
    printf("%s : Compression ratio: %lf with error bound %lf\n", compressor_string.c_str(), cr, eb);
    return cr;
}

double get_real_cr(const std::string _, pressio_data *input, const double eb) {
    if (is_reversed_cmp_dimension[_]) {
        pressio_data_reshape(input, reversed_dims.size(), reversed_dims.data());
    }
    else {
        pressio_data_reshape(input, dims.size(), dims.data());
    }

    std::string compressor_string = _;
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    }

    pressio_options compression_options{
        {"sz3:metric", "composite"},
        {"composite:plugins", std::vector{"size"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

    if(compressor->compress(input, &compressed) != 0) {
        std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    size_t compressed_size = compressed.size_in_bytes();
    return compressed_size;
}

size_t compress_to_destination(const std::string _, pressio_data *input, const double eb, char *&output) {

    if (is_reversed_cmp_dimension[_]) {
        pressio_data_reshape(input, reversed_dims.size(), reversed_dims.data());
    }
    else {
        pressio_data_reshape(input, dims.size(), dims.data());
    }

    std::string compressor_string = _;
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    }

    pressio_options compression_options{
        {"sz3:metric", "composite"},
        {"composite:plugins", std::vector{"size"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

    if(compressor->compress(input, &compressed) != 0) {
        std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    size_t compressed_size = compressed.size_in_bytes();
    output = reinterpret_cast<char *>(compressed.release());
    output = (char*)std::realloc(output, compressed_size + sizeof(ALGO));
    memcpy(output + compressed_size, &compressor_string_to_enum[compressor_string], sizeof(ALGO));
    return compressed_size;
}

#endif