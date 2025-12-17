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
#include <limits>

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

#include "Timer.hpp"

pressio library;

using namespace std::string_literals;

static const auto compressor_list = {"sz3", "zfp", "szx", "sperr"};
// static const auto compressor_list = {"sz3", "zfp", "szx"};
// static const auto compressor_list = {"sperr"};

std::map<std::string, uint8_t> is_reversed_surrogate_dimension = {
    {"sz3", 0}, 
    {"zfp", 1}, 
    {"szx", 0}, 
    {"sperr", 0}, 
    {"szp", 0}
};

std::map<std::string, ALGO> compressor_string_to_enum = {
    {"sz3", ALGO::ALGO_SZ3}, 
    {"zfp", ALGO::ALGO_ZFP}, 
    {"szx", ALGO::ALGO_SZx}, 
    {"sperr", ALGO::ALGO_SPERR}, 
    {"szp", ALGO::ALGO_SZp}
};

std::map<ALGO, std::string> compressor_enum_to_string = {
    {ALGO::ALGO_SZ3, "sz3"}, 
    {ALGO::ALGO_ZFP, "zfp"}, 
    {ALGO::ALGO_SZx, "szx"}, 
    {ALGO::ALGO_SPERR, "sperr"}, 
    {ALGO::ALGO_SZp, "szp"}
};

std::map<std::string, uint8_t> &is_reversed_dimension = is_reversed_surrogate_dimension;

std::map<std::string, uint8_t> is_reversed_cmp_dimension = {
    {"sz3", 0}, 
    {"zfp", 0}, 
    {"szx", 0}, 
    {"sperr", 0}, 
    {"szp", 0}
};

std::map<std::string, std::string> surrogate_option = {
    {"sz3", "sz3_surrogate:abs_error_bound"}, 
    {"zfp", "zfp_surrogate:accuracy"}, 
    {"szx", "szx_surrogate:abs_err_bound"}, 
    {"sperr", "sperr_surrogate:abs_err_bound"}, 
    {"szp", "szp_surrogate:abs_err_bound"}
};

double get_surrogate_cr(const std::string _, pressio_data *input, const double eb);
double get_surrogate_psnr(const std::string _, pressio_data *input, const double eb);
double get_surrogate_ssim(const std::string _, pressio_data *input, const double eb);

std::pair<double, double> get_surrogate_cr_with_psnr(const std::string _, pressio_data *input, const double psnr);
double get_surrogate_eb_with_cr(const std::string _, pressio_data *input, const double cr);
std::pair<double, double> get_surrogate_psnr_with_cr(const std::string _, pressio_data *input, const double cr);

double get_real_cr(const std::string _, pressio_data *input, const double eb);
double get_real_psnr(const std::string _, pressio_data *input, const double eb);
double get_real_ssim(const std::string _, pressio_data *input, const double eb);
size_t compress_to_destination(const std::string _, pressio_data *input, const double eb, char *&output);

size_t compress_fix_eb_best_cr(pressio_data *input, const double eb, char *&output);
size_t compress_fix_eb_best_psnr(pressio_data *input, const double eb, char *&output);
// size_t compress_fix_eb_best_ssim(pressio_data *input, const double eb, char *&output);

size_t compress_fix_psnr_best_cr(pressio_data *input, const double psnr, char *&output);
// size_t compress_fix_ssim_best_cr(pressio_data *input, const double ssim, char *&output);

size_t compress_fix_cr_best_eb(pressio_data *input, const double cr, char *&output);
size_t compress_fix_cr_best_psnr(pressio_data *input, const double cr, char *&output);
// size_t compress_fix_cr_best_ssim(pressio_data *input, const double cr, char *&output);

using function_t = std::function<size_t(
    pressio_data *, 
    const double, 
    char *&)>;

std::unordered_map<std::string, function_t> func_mp = {
    {"eb_cr", compress_fix_eb_best_cr}, 
    {"eb_psnr", compress_fix_eb_best_psnr}, 
    {"psnr_cr", compress_fix_psnr_best_cr}, 
    // {"ssim_cr", compress_fix_ssim_best_cr}, 
    {"cr_eb", compress_fix_cr_best_eb}, 
    {"cr_psnr", compress_fix_cr_best_psnr}, 
    // {"cr_ssim", compress_fix_cr_best_ssim}
};

std::string metrics[4] = {"eb", "psnr", "ssim", "cr"};

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

    std::string method_string = metrics[conf.method_id / 4] + "_" + metrics[conf.method_id % 4];

    return func_mp[method_string](&input, compression_parameter, compressedData);
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

size_t compress_fix_eb_best_psnr(pressio_data *input, const double eb, char *&output) {
    
    double best_psnr = 0;
    std::string best_compressor = "null";
    // test the best compressor
    for (const auto &_ : compressor_list) {
        auto current_psnr = get_surrogate_psnr(_, input, eb);
        if (current_psnr > best_psnr) {
            best_psnr = current_psnr;
            best_compressor = _;
        }
    }

    // std::cout << "Best compressor: " << best_compressor << " with CR: " << best_cr << std::endl;

    size_t compressed_size = compress_to_destination(best_compressor, input, eb, output);
    return compressed_size;
}

size_t compress_fix_eb_best_ssim(pressio_data *input, const double eb, char *&output) {
    
    double best_ssim = 0;
    std::string best_compressor = "null";
    // test the best compressor
    for (const auto &_ : compressor_list) {
        auto current_ssim = get_surrogate_ssim(_, input, eb);
        auto real_ssim = get_real_ssim(_, input, eb);
        if (current_ssim > best_ssim) {
            best_ssim = current_ssim;
            best_compressor = _;
        }
    }

    // std::cout << "Best compressor: " << best_compressor << " with CR: " << best_cr << std::endl;

    size_t compressed_size = compress_to_destination(best_compressor, input, eb, output);
    return compressed_size;
}

size_t compress_fix_psnr_best_cr(pressio_data *input, const double psnr, char *&output) {

    double best_eb = std::numeric_limits<double>::infinity();
    double best_cr = 0;
    std::string best_compressor = "null";

    for (const auto &_ : compressor_list) {
        auto result = get_surrogate_cr_with_psnr(_, input, psnr);
        auto &current_eb = result.first;
        auto &current_cr = result.second;
        if (current_cr > best_cr) {
            best_eb = current_eb;
            best_cr = current_cr;
            best_compressor = _;
        }
    }

    // std::cout << "best_compressor = " << best_compressor << " eb = " << best_eb << std::endl;

    size_t compressed_size = compress_to_destination(best_compressor, input, best_eb, output);
    return compressed_size;
}


size_t compress_fix_cr_best_eb(pressio_data *input, const double cr, char *&output) {
    double best_eb = std::numeric_limits<double>::infinity();
    std::string best_compressor = "null";

    for (const auto &_ : compressor_list) {
        double current_eb = get_surrogate_eb_with_cr(_, input, cr);
        if (current_eb < best_eb) {
            best_eb = current_eb;
            best_compressor = _;
        }
    }

    // std::cout << "best_compressor = " << best_compressor << " eb = " << best_eb << std::endl;

    size_t compressed_size = compress_to_destination(best_compressor, input, best_eb, output);
    return compressed_size;
}

size_t compress_fix_cr_best_psnr(pressio_data *input, const double cr, char *&output) {

    double best_eb = std::numeric_limits<double>::infinity();
    double best_psnr = 0;
    std::string best_compressor = "null";

    for (const auto &_ : compressor_list) {
        auto result = get_surrogate_psnr_with_cr(_, input, cr);
        auto &current_eb = result.first;
        auto &current_psnr = result.second;
        if (current_psnr > best_psnr) {
            best_eb = current_eb;
            best_psnr = current_psnr;
            best_compressor = _;
        }
        // printf("compressor = %s, eb = %lf, cr = %lf, psnr = %lf\n", _, current_eb, get_surrogate_cr(_, input, current_eb), current_psnr);
    }

    // std::cout << "best_compressor = " << best_compressor << " eb = " << best_eb << std::endl;

    size_t compressed_size = compress_to_destination(best_compressor, input, best_eb, output);
    return compressed_size;
}

double get_range(pressio_data *input) {
    size_t n = input -> num_elements();
    void *data = input -> data();
    if (input -> dtype() == pressio_float_dtype) {
        float *vec = reinterpret_cast<float *>(data);
        float __max = vec[0];
        for (size_t i = 1; i < n; i++) __max = std::max(__max, vec[i]);
        return static_cast<double>(__max);
    }
    else if (input -> dtype() == pressio_double_dtype) {
        double *vec = reinterpret_cast<double *>(data);
        double __max = vec[0];
        for (size_t i = 1; i < n; i++) __max = std::max(__max, vec[i]);
        return __max;
    }
    else {
        exit(-1);
    }
    return 0;
}

std::pair<double, double> get_surrogate_cr_with_psnr(const std::string _, pressio_data *input, const double target_psnr) {

    double l = -16, r = log10(get_range(input)), mid;
    while((r - l) > 1e-1) {
        mid = (l + r) / 2.;
        double current_psnr = get_surrogate_psnr(_, input, pow(10, mid));
        if (current_psnr > target_psnr) {
            l = mid;
        } else {
            r = mid;
        }
        // printf("eb = %lf, psnr = %lf\n", pow(10, mid), current_psnr);
    }

    // printf("real psnr = %lf\n", get_real_psnr(_, input, pow(10, mid)));

    return {pow(10, mid), get_surrogate_cr(_, input, pow(10, mid))};
}

double get_surrogate_eb_with_cr(const std::string _, pressio_data *input, const double target_cr) {

    double l = -16, r = log10(get_range(input)), mid;
    while((r - l) > 1e-1) {
        mid = (l + r) / 2.;
        double current_cr = get_surrogate_cr(_, input, pow(10, mid));
        if (current_cr < target_cr) {
            l = mid;
        } else {
            r = mid;
        }
        // printf("eb = %lf, cr = %lf\n", pow(10, mid), current_cr);
    }

    return pow(10, mid);
}

std::pair<double, double> get_surrogate_psnr_with_cr(const std::string _, pressio_data *input, const double target_cr) {
    
    double l = -16, r = log10(get_range(input)), mid;
    while((r - l) > 1e-1) {
        mid = (l + r) / 2.;
        double current_cr = get_surrogate_cr(_, input, pow(10, mid));
        if (current_cr < target_cr) {
            l = mid;
        } else {
            r = mid;
        }
        // printf("eb = %lf, cr = %lf\n", pow(10, mid), current_cr);
    }

    return {pow(10, mid), get_surrogate_psnr(_, input, pow(10, mid))};
}

double get_surrogate_cr(const std::string _, pressio_data *input, const double eb) {

    SZ3::Timer timer(true);

    if (is_reversed_surrogate_dimension[_]) {
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
        // {_ + "_surrogate:abs_error_bound", eb}, 
        // {_ + "_surrogate:error_bound_mode", "ABS"}, 
        {surrogate_option[_], eb}, 
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
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double time = timer.stop();
    // printf("Surrogate time: %lf seconds\n", time);

    double cr;
    assert_defined((_ + "_surrogate:cr").c_str(), &cr);
    // printf("%s : Compression ratio: %lf with error bound %lf\n", compressor_string.c_str(), cr, eb);
    return cr;
}

double get_real_cr(const std::string _, pressio_data *input, const double eb) {

    SZ3::Timer timer(true);

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
        {compressor_string + ":metric", "composite"},
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

    auto metrics = compressor->get_metrics_results();
    auto assert_defined = [&](const char* key, auto value){
        if(metrics.get(key, value)!= pressio_options_key_set) {
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double time = timer.stop();
    // printf("Compression time: %lf seconds\n", time);

    double cr;
    assert_defined("size:compression_ratio", &cr);
    // printf("%s : Compression ratio: %lf with error bound %lf\n", compressor_string.c_str(), cr, eb);

    return cr;
}

double get_surrogate_psnr(const std::string _, pressio_data *input, const double eb) {

    SZ3::Timer timer(true);

    if (is_reversed_surrogate_dimension[_]) {
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
        // {_ + "_surrogate:abs_error_bound", eb}, 
        // {_ + "_surrogate:error_bound_mode", "ABS"}, 
        {surrogate_option[_], eb}, 
        {_ + "_surrogate:metric", "psnr"}
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
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double time = timer.stop();
    // printf("Surrogate time: %lf seconds\n", time);

    double psnr;
    assert_defined((_ + "_surrogate:psnr").c_str(), &psnr);
    // printf("%s : estimate PSNR: %lf with error bound %lf\n", compressor_string.c_str(), psnr, eb);
    return psnr;
}

double get_real_psnr(const std::string _, pressio_data *input, const double eb) {

    SZ3::Timer timer(true);

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
        {compressor_string + ":metric", "composite"},
        {"composite:plugins", std::vector{"error_stat"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});
    pressio_data output = pressio_data::owning(pressio_data_dtype(input), dims);

    if(compressor->compress(input, &compressed) != 0) {
        std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    if(compressor->decompress(&compressed, &output) != 0) {
        std::cerr << compressor->error_msg() << std::endl;
        exit(compressor->error_code());
    }

    auto metrics = compressor->get_metrics_results();
    auto assert_defined = [&](const char* key, auto value){
        if(metrics.get(key, value)!= pressio_options_key_set) {
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double time = timer.stop();
    // printf("Compression time: %lf seconds\n", time);

    double psnr;
    assert_defined("error_stat:psnr", &psnr);
    // printf("%s : PSNR: %lf with error bound %lf\n", compressor_string.c_str(), psnr, eb);

    return psnr;
}

double get_surrogate_ssim(const std::string _, pressio_data *input, const double eb) {

    SZ3::Timer timer(true);

    if (is_reversed_surrogate_dimension[_]) {
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
        {surrogate_option[_], eb}, 
        {_ + "_surrogate:metric", "ssim"}
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
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double time = timer.stop();
    printf("Surrogate time: %lf seconds\n", time);

    double ssim;
    assert_defined((_ + "_surrogate:ssim").c_str(), &ssim);
    printf("%s : estimate SSIM: %lf with error bound %lf\n", compressor_string.c_str(), ssim, eb);
    return ssim;
}

double get_real_ssim(const std::string _, pressio_data *input, const double eb) {

    SZ3::Timer timer(true);

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
        {compressor_string + ":metric", "composite"},
        {"composite:plugins", std::vector{"ssim"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});
    pressio_data output = pressio_data::owning(pressio_data_dtype(input), dims);

    if(compressor->compress(input, &compressed) != 0) {
        std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    if(compressor->decompress(&compressed, &output) != 0) {
        std::cerr << compressor->error_msg() << std::endl;
        exit(compressor->error_code());
    }

    auto metrics = compressor->get_metrics_results();
    auto assert_defined = [&](const char* key, auto value){
        if(metrics.get(key, value)!= pressio_options_key_set) {
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double time = timer.stop();
    printf("Compression time: %lf seconds\n", time);

    double ssim;
    assert_defined("ssim:ssim", &ssim);
    printf("%s : SSIM: %lf with error bound %lf\n", compressor_string.c_str(), ssim, eb);

    return ssim;
}

size_t compress_to_destination(const std::string _, pressio_data *input, const double eb, char *&output) {

    SZ3::Timer timer(true);

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
        {compressor_string + ":metric", "composite"},
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

    output = reinterpret_cast<char *>(std::malloc(compressed_size + sizeof(ALGO) + sizeof(size_t) + sizeof(double)));
    char *tail = output;
    memcpy(tail, &compressed_size, sizeof(size_t));
    tail += sizeof(size_t);
    memcpy(tail, &compressor_string_to_enum[compressor_string], sizeof(ALGO));
    tail += sizeof(ALGO);
    memcpy(tail, &eb, sizeof(double));
    tail += sizeof(double);
    memcpy(tail, compressed.data(), compressed_size);

    double time = timer.stop();
    // printf("Compression time: %lf seconds\n", time);

    return compressed_size;
}

template <typename T>
void ROCCI_select_decompress(Config &conf, T *decompressed_data, char *compressed_data) {

    size_t compressed_size;
    memcpy(&compressed_size, compressed_data, sizeof(size_t));
    compressed_data += sizeof(size_t);
    ALGO compressor_enum;
    memcpy(&compressor_enum, compressed_data, sizeof(ALGO));
    compressed_data += sizeof(ALGO);
    std::string compressor_string = compressor_enum_to_string[compressor_enum];

    double eb;
    memcpy(&eb, compressed_data, sizeof(double));
    compressed_data += sizeof(double);

    constexpr auto float_type = std::is_same<T, float>::value ? pressio_float_dtype : pressio_double_dtype;
    pressio_data input = pressio_data::nonowning(pressio_byte_dtype, compressed_data, {compressed_size});
    pressio_data output = pressio_data::owning(float_type, conf.dims);

    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    }

    pressio_options compression_options{
        {"pressio:metric", "composite"s},
        {"composite:plugins", std::vector{"size"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    if (compressor->decompress(&input, &output) != 0) {
        std::cerr << "Decompression error: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    memcpy(decompressed_data, output.data(), conf.num * sizeof(T));

    return;
}

#endif