#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>

#include <libpressio.h>
#include <libpressio_ext/cpp/libpressio.h>
#include <libpressio_ext/cpp/compressor.h>
#include <libpressio_ext/cpp/options.h>
#include <libpressio_ext/cpp/data.h>
#include <libpressio_meta.h>

pressio library;

using namespace std::string_literals;

double eb = 1e-2;

void testSZ3_surrogate(pressio_data *input);
void testSZ3(pressio_data *input);

void testZFP_surrogate(pressio_data *input);
void testZFP(pressio_data *input);

void testSZx_surrogate(pressio_data *input);
void testSZx(pressio_data *input);


signed main(int argc, char* argv[]) {

    libpressio_register_all();

    // std::cout << library.supported_compressors() << std::endl;
    // return 0;

    // For data and comp info
    std::string input_file = "/anvil/projects/x-cis240669/kai/data/hurricane-100x500x500/Uf48.bin.dat";
    pressio_dtype type = pressio_float_dtype;

    // std::vector<size_t> dims {100, 500, 500};
    std::vector<size_t> dims {500, 500, 100};
    
    // Read the input data

    pressio_data metadata = pressio_data::owning(type, dims);

    pressio_io io = library.get_io("posix");
    if(!io) {
        std::cerr << library.err_msg() << std::endl;
        exit(library.err_code());
    }
    if(io->set_options({
      {"io:path", input_file},
    })) {
        std::cerr << io->error_msg() << std::endl;
        exit(io->error_code());
    }
    pressio_data* input = io->read(&metadata);
    // pressio_data* input_reverse_dims = 
    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});
    pressio_data output = pressio_data::owning(type, dims);

    // testSZ3_surrogate(input);
    // testSZ3(input);

    testZFP_surrogate(input);
    testZFP(input);

    // testSZ3_surrogate(input);
    // testSZ3(input);
    













}

std::string crs =  "cr";

void testSZ3_surrogate(pressio_data *input) {

    std::string compressor_string = "sz3_surrogate";
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
    }

    pressio_options compression_options{
        {"sz3_surrogate:abs_error_bound", eb}, 
        {"sz3_surrogate:metric", "cr"}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
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
    assert_defined(("sz3_surrogate:" + crs).c_str(), &cr);
    printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
}

void testSZ3(pressio_data *input) {

    std::string compressor_string = "sz3";
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
    }

    pressio_options compression_options{
        {"sz3:metric", "composite"},
        {"composite:plugins", std::vector{"size"s, "error_stat"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
    }

    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

    if(compressor->compress(input, &compressed) != 0) {
        std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
        exit(-1);
    }

    // pressio_data output = pressio_data::owning(pressio_float_dtype, {500, 500, 100});

    // if(compressor->decompress(&compressed, &output) != 0) {
    //     std::cerr << compressor->error_msg() << std::endl;
    //     exit(compressor->error_code());
    // }

    auto metrics = compressor->get_metrics_results();
    auto assert_defined = [&](const char* key, auto value){
        if(metrics.get(key, value)!= pressio_options_key_set) {
            throw std::runtime_error(key + " was note set"s);
        }
    };

    double cr;
    // assert_defined("error_stat:psnr", &cr);
    assert_defined("size:compression_ratio", &cr);
    printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
}

void testZFP_surrogate(pressio_data *input) {

    std::string compressor_string = "zfp_surrogate";
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
    }

    pressio_options compression_options{
        {"zfp_surrogate:abs_error_bound", eb}, 
        {"zfp_surrogate:metric", crs}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
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
    assert_defined(("zfp_surrogate:" + crs).c_str(), &cr);
    printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
}

void testZFP(pressio_data *input) {

    std::string compressor_string = "zfp";
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
    }

    pressio_options compression_options{
        {"zfp:metric", "composite"},
        {"composite:plugins", std::vector{"size"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
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
    assert_defined("size:compression_ratio", &cr);
    printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
}

void testSZx_surrogate(pressio_data *input) {

    std::string compressor_string = "szx_surrogate";
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
    }

    pressio_options compression_options{
        {"szx_surrogate:abs_error_bound", eb}, 
        {"szx_surrogate:metric", crs}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
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
    assert_defined(("szx_surrogate:" + crs).c_str(), &cr);
    printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
}

void testSZx(pressio_data *input) {

    std::string compressor_string = "szx";
    auto compressor = library.get_compressor(compressor_string);
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
    }

    pressio_options compression_options{
        {"szx:metric", "composite"},
        {"composite:plugins", std::vector{"size"s}},
        {"pressio:abs", eb}
    };

    if(compressor->set_options(compression_options) != 0) {
        std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
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
    assert_defined("size:compression_ratio", &cr);
    printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
}