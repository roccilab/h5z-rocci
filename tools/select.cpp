#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>

#include <libpressio.h>
#include <libpressio_ext/cpp/libpressio.h>
#include <libpressio_ext/cpp/compressor.h>
#include <libpressio_ext/cpp/options.h>
#include <libpressio_ext/cpp/data.h>
#include <libpressio_meta.h>

pressio library;

using namespace std::string_literals;

static const auto compressor_list = {"sz3", "zfp", "szx"};

double get_surrogate_cr(const std::string _, const pressio_data *input, const double eb);
void compress_to_destination(const std::string _, const pressio_data *input, const double eb, const char *output_file);

void compress_fix_eb_best_cr(const pressio_data *input, const double eb, const char *output_file);
void compress_fix_psnr_best_cr(const pressio_data *input, const double psnr, const char *output_file);
void compress_fix_ssim_best_cr(const pressio_data *input, const double ssim, const char *output_file);

void compress_fix_cr_best_eb(const pressio_data *input, const double cr, const char *output_file);
void compress_fix_cr_best_psnr(const pressio_data *input, const double cr, const char *output_file);
void compress_fix_cr_best_ssim(const pressio_data *input, const double cr, const char *output_file);

/*
./tools/select \
/anvil/projects/x-cis240669/kai/data/hurricane-100x500x500/Uf48.bin.dat \
3 500 500 100 \
~/scratch/tem eb_cr \
1e-1
*/

using function_t = std::function<void(
    const pressio_data *, 
    const double, 
    const char *)>;

std::unordered_map<std::string, function_t> func_mp = {
    {"eb_cr", compress_fix_eb_best_cr}, 
    {"psnr_cr", compress_fix_psnr_best_cr}, 
    {"ssim_cr", compress_fix_ssim_best_cr}, 
    {"cr_eb", compress_fix_cr_best_eb}, 
    {"cr_psnr", compress_fix_cr_best_psnr}, 
    {"cr_ssim", compress_fix_cr_best_ssim}
};

signed main(int argc, char* argv[]) {

    char *input_file, *output_file;
    std::vector<size_t> dims;
    double compression_parameter; //eb, psnr, ssim or cr, which are all in double type
    pressio_dtype type = pressio_float_dtype;

    function_t function;

    {
        int i = 0;
        input_file = argv[++i];
        dims.resize(std::stoi(argv[++i]));
        for (auto &dim : dims) {
            dim = std::stoll(argv[++i]);
        }
        output_file = argv[++i];
        std::string mode = std::string(argv[++i]);
        if (func_mp.find(mode) == func_mp.end()) {
            std::cerr << "Invalid mode: " << mode << std::endl;
            std::cerr << "Valid modes are: ";
            for (const auto &func : func_mp) {
                std::cerr << func.first << " ";
            }
            std::cerr << std::endl;
            exit(-1);
        }
        function = func_mp[mode];
        compression_parameter = std::stod(argv[++i]);
    }

    libpressio_register_all();
    
    // Read the input data

    pressio_data metadata = pressio_data::owning(type, dims);

    pressio_io io = library.get_io("posix");
    if(!io) {
        std::cerr << library.err_msg() << std::endl;
        exit(library.err_code());
    }
    if(io->set_options({
      {"io:path", std::string(input_file)},
    })) {
        std::cerr << io->error_msg() << std::endl;
        exit(io->error_code());
    }
    pressio_data* input = io->read(&metadata);
    pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});
    pressio_data output = pressio_data::owning(type, dims);

    // testSZx_surrogate(input, eb);
    // testSZx(input, eb);

    function(input, compression_parameter, output_file);
}

void compress_fix_eb_best_cr(const pressio_data *input, const double eb, const char *output_file) {

    double best_cr = 0;
    std::string best_compressor = "null";
    for (const auto &_ : compressor_list) {
        auto current_cr = get_surrogate_cr(_, input, eb);
        if (current_cr > best_cr) {
            best_cr = current_cr;
            best_compressor = _;
        }
    }

    std::cout << "Best compressor: " << best_compressor << " with CR: " << best_cr << std::endl;
    
    compress_to_destination(best_compressor, input, eb, output_file);
}

void compress_fix_psnr_best_cr(const pressio_data *input, const double psnr, const char *output_file) {


}

void compress_fix_ssim_best_cr(const pressio_data *input, const double ssim, const char *output_file) {

}

void compress_fix_cr_best_eb(const pressio_data *input, const double cr, const char *output_file) {

    double best_eb = -1;
    std::string best_compressor = "null";

    for (const auto &_ : compressor_list) {

        double l = 0, r = 5, mid;
        // 10 ^ -

        while(r - l > 1e-2) {
            mid = (l + r) / 2.;
            auto current_cr = get_surrogate_cr(_, input, pow(10, -mid));
            if (current_cr > cr) {
                l = mid;
            } else {
                r = mid;
            }
        }

        if (pow(10, -l) > best_eb) {
            best_eb = pow(10, -l);
            best_compressor = _;
        }

        std::cout << "Compressor: " << _ << " with EB: " << pow(10, -l) << std::endl;
    }

    std::cout << "Best compressor: " << best_compressor << " with EB: " << best_eb << std::endl;
}

void compress_fix_cr_best_psnr(const pressio_data *input, const double cr, const char *output_file) {

}

void compress_fix_cr_best_ssim(const pressio_data *input, const double cr, const char *output_file) {

}


double get_surrogate_cr(const std::string _, const pressio_data *input, const double eb) {

    const std::string &compressor_string = _;
    auto compressor = library.get_compressor(_ + "_surrogate");
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
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

    /*
    tem
    */

    /*
    tem
    */

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

double get_surrogate_psnr(const std::string _, const pressio_data *input, const double eb) {

    const std::string &compressor_string = _;
    auto compressor = library.get_compressor(_ + "_surrogate");
    if (!compressor) {
        std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
        exit(-1);
    }

    pressio_options compression_options{
        {_ + "_surrogate:abs_error_bound", eb}, 
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
            throw std::runtime_error(key + " was note set"s);
        }
    };

    double psnr;
    assert_defined((_ + "_surrogate:psnr").c_str(), &psnr);
    // printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
    return psnr;
}

void compress_to_destination(const std::string _, const pressio_data *input, const double eb, const char *output_file) {

    std::string compressor_string = "sz3";
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

    auto metrics = compressor->get_metrics_results();
    auto assert_defined = [&](const char* key, auto value){
        if(metrics.get(key, value)!= pressio_options_key_set) {
            throw std::runtime_error(key + " was note set"s);
        }
    };

    double cr;
    assert_defined("size:compression_ratio", &cr);
    printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);


    pressio_io io = library.get_io("posix");
    if(!io) {
        std::cerr << library.err_msg() << std::endl;
        exit(library.err_code());
    }
    if(io->set_options({
      {"io:path", std::string(output_file)},
    })) {
        std::cerr << io->error_msg() << std::endl;
        exit(io->error_code());
    }

    if(io->write(&compressed) != 0) {
        std::cerr << "Write error: " << io->error_msg() << std::endl;
        exit(-1);
    } else {
        std::cerr << "Successfully wrote compressed data to " << output_file << std::endl;
    }

    
}


// old test functions

// void testSZ3_surrogate(pressio_data *input, double eb) {

//     std::string compressor_string = "sz3_surrogate";
//     auto compressor = library.get_compressor(compressor_string);
//     if (!compressor) {
//         std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
//     }

//     pressio_options compression_options{
//         {"sz3_surrogate:abs_error_bound", eb}, 
//         {"sz3_surrogate:metric", "cr"}
//     };

//     if(compressor->set_options(compression_options) != 0) {
//         std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
//     }

//     pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

//     if(compressor->compress(input, &compressed) != 0) {
//         std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     }

//     auto metrics = compressor->get_metrics_results();
//     auto assert_defined = [&](const char* key, auto value){
//         if(metrics.get(key, value)!= pressio_options_key_set) {
//             throw std::runtime_error(key + " was note set"s);
//         }
//     };

//     double cr;
//     assert_defined("sz3_surrogate:cr", &cr);
//     printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
// }

// void testSZ3(pressio_data *input, double eb) {

//     std::string compressor_string = "sz3";
//     auto compressor = library.get_compressor(compressor_string);
//     if (!compressor) {
//         std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
//     }

//     pressio_options compression_options{
//         {"sz3:metric", "composite"},
//         {"composite:plugins", std::vector{"size"s}},
//         {"pressio:abs", eb}
//     };

//     if(compressor->set_options(compression_options) != 0) {
//         std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
//     }

//     pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

//     if(compressor->compress(input, &compressed) != 0) {
//         std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     }

//     auto metrics = compressor->get_metrics_results();
//     auto assert_defined = [&](const char* key, auto value){
//         if(metrics.get(key, value)!= pressio_options_key_set) {
//             throw std::runtime_error(key + " was note set"s);
//         }
//     };

//     double cr;
//     assert_defined("size:compression_ratio", &cr);
//     printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
// }

// void testZFP_surrogate(pressio_data *input, double eb) {

//     std::string compressor_string = "zfp_surrogate";
//     auto compressor = library.get_compressor(compressor_string);
//     if (!compressor) {
//         std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
//     }

//     pressio_options compression_options{
//         {"zfp_surrogate:abs_error_bound", eb}, 
//         {"zfp_surrogate:metric", "cr"}
//     };

//     if(compressor->set_options(compression_options) != 0) {
//         std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
//     }

//     pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

//     if(compressor->compress(input, &compressed) != 0) {
//         std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     }

//     auto metrics = compressor->get_metrics_results();
//     auto assert_defined = [&](const char* key, auto value){
//         if(metrics.get(key, value)!= pressio_options_key_set) {
//             throw std::runtime_error(key + " was note set"s);
//         }
//     };

//     double cr;
//     assert_defined("zfp_surrogate:cr", &cr);
//     printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
// }

// void testZFP(pressio_data *input, double eb) {

//     std::string compressor_string = "zfp";
//     auto compressor = library.get_compressor(compressor_string);
//     if (!compressor) {
//         std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
//     }

//     pressio_options compression_options{
//         {"zfp:metric", "composite"},
//         {"composite:plugins", std::vector{"size"s}},
//         {"pressio:abs", eb}
//     };

//     if(compressor->set_options(compression_options) != 0) {
//         std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
//     }

//     pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

//     if(compressor->compress(input, &compressed) != 0) {
//         std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     }

//     auto metrics = compressor->get_metrics_results();
//     auto assert_defined = [&](const char* key, auto value){
//         if(metrics.get(key, value)!= pressio_options_key_set) {
//             throw std::runtime_error(key + " was note set"s);
//         }
//     };

//     double cr;
//     assert_defined("size:compression_ratio", &cr);
//     printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
// }

// void testSZx_surrogate(pressio_data *input, double eb) {

//     std::string compressor_string = "szx_surrogate";
//     auto compressor = library.get_compressor(compressor_string);
//     if (!compressor) {
//         std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
//     }

//     pressio_options compression_options{
//         {"szx_surrogate:abs_error_bound", eb}, 
//         {"szx_surrogate:metric", "cr"}
//     };

//     if(compressor->set_options(compression_options) != 0) {
//         std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
//     }

//     pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

//     if(compressor->compress(input, &compressed) != 0) {
//         std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     }

//     auto metrics = compressor->get_metrics_results();
//     auto assert_defined = [&](const char* key, auto value){
//         if(metrics.get(key, value)!= pressio_options_key_set) {
//             throw std::runtime_error(key + " was note set"s);
//         }
//     };

//     double cr;
//     assert_defined("szx_surrogate:cr", &cr);
//     printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
// }

// void testSZx(pressio_data *input, double eb) {

//     std::string compressor_string = "szx";
//     auto compressor = library.get_compressor(compressor_string);
//     if (!compressor) {
//         std::cerr << "Failed to load " << compressor_string << " compressor" << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully load " << compressor_string << " compressor" << std::endl;
//     }

//     pressio_options compression_options{
//         {"szx:metric", "composite"},
//         {"composite:plugins", std::vector{"size"s}},
//         {"pressio:abs", eb}
//     };

//     if(compressor->set_options(compression_options) != 0) {
//         std::cerr << "Failed to set options: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     } else {
//         std::cerr << "Successfully set options: " << compressor->error_msg() << std::endl;
//     }

//     pressio_data compressed = pressio_data::empty(pressio_byte_dtype, {});

//     if(compressor->compress(input, &compressed) != 0) {
//         std::cerr << "Compression error: " << compressor->error_msg() << std::endl;
//         exit(-1);
//     }

//     auto metrics = compressor->get_metrics_results();
//     auto assert_defined = [&](const char* key, auto value){
//         if(metrics.get(key, value)!= pressio_options_key_set) {
//             throw std::runtime_error(key + " was note set"s);
//         }
//     };

//     double cr;
//     assert_defined("size:compression_ratio", &cr);
//     printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
// }