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

// static const auto compressor_list = {"sz3", "zfp", "szx"};
static const auto compressor_list = {"szx"};

double get_surrogate_cr(const std::string _, pressio_data *input, const double eb);
double get_real_cr(const std::string _, pressio_data *input, const double eb);
void compress_to_destination(const std::string _, pressio_data *input, const double eb, const char *output_file);

void compress_fix_eb_best_cr(pressio_data *input, const double eb, const char *output_file);
void compress_fix_psnr_best_cr(pressio_data *input, const double psnr, const char *output_file);
void compress_fix_ssim_best_cr(pressio_data *input, const double ssim, const char *output_file);

void compress_fix_cr_best_eb(pressio_data *input, const double cr, const char *output_file);
void compress_fix_cr_best_psnr(pressio_data *input, const double cr, const char *output_file);
void compress_fix_cr_best_ssim(pressio_data *input, const double cr, const char *output_file);

void print_list(auto &list) {
    std::cout << "{";
    for (int i = 0; i < list.size(); i++) {
        std::cout << list[i];
        if (i != list.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "}" << std::endl;
}

void print_list(auto &eb_list, auto &cr_list) {
    print_list(eb_list);
    print_list(cr_list);
}

void generate_compression_ratio_list(pressio_data *input, std::vector<double> &eb_list) {
    std::vector<double> cr_list;
    for (auto &eb : eb_list) {
        cr_list.push_back(get_surrogate_cr("sz3", input, eb));
    }
    print_list(eb_list, cr_list);

    cr_list.clear();
    for (auto &eb : eb_list) {
        cr_list.push_back(get_real_cr("sz3", input, eb));
    }
    print_list(eb_list, cr_list);

}

void generate_compression_ratio_list_linear(pressio_data *input) {
    std::vector<double> eb_list;
    for (int i = 0; i <= 15; i++) {
        eb_list.push_back(1. * (i + 1) / 16);
    }
    generate_compression_ratio_list(input, eb_list);
}

void generate_compression_ratio_list_log(pressio_data *input) {
    std::vector<double> eb_list;
    for (int i = 0; i <= 15; i++) {
        eb_list.push_back(pow(10, -1.0 * i / 3));
    }
    generate_compression_ratio_list(input, eb_list);
}

void generate_fix_compression_ratio_error_list(pressio_data *input) {
    std::vector<double> cr_list;
    // for (int i = 2; i <= 20; i++) cr_list.push_back(i);
    for (int i = 12; i <= 20; i++) cr_list.push_back(i / 10.0);
    std::vector<double> cr_surrogate_list, cr_real_list;

    std::string _ = "zfp";

    for (auto &cr : cr_list) {
        double l = -1, r = 5, mid;

        while(r - l > 1e-2) {
            mid = (l + r) / 2.;
            printf("%lf\n", pow(10, -mid));
            auto current_cr = get_surrogate_cr(_, input, pow(10, -mid));
            if (current_cr > cr) {
                l = mid;
            } else {
                r = mid;
            }
        }

        cr_surrogate_list.push_back(get_surrogate_cr(_, input, pow(10, -l)));
        cr_real_list.push_back(get_real_cr(_, input, pow(10, -l)));
    }

    print_list(cr_list);
    print_list(cr_surrogate_list, cr_real_list);
}

/*
./tools/select \
/anvil/projects/x-cis240669/kai/data/hurricane-100x500x500/Uf48.bin.dat \
3 500 500 100 \
~/scratch/tem eb_cr \
1e-1

./tools/select \
/anvil/projects/x-cis240669/kai/data/hurricane-100x500x500/Uf48.bin.dat \
3 500 500 100 \
~/scratch/tem cr_eb \
10
*/

std::map<std::string, uint8_t> is_reversed_surrogate_dimension = {
    {"sz3", 0}, 
    {"zfp", 0}, 
    {"szx", 0}, 
    {"sperr", 0}, 
    {"szp", 1}
};
std::map<std::string, uint8_t> &is_reversed_dimension = is_reversed_surrogate_dimension;

std::vector<size_t> dims, reversed_dims;

using function_t = std::function<void(
    pressio_data *, 
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

    reversed_dims = dims;
    std::reverse(reversed_dims.begin(), reversed_dims.end());

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

    // function(input, compression_parameter, output_file);
    // generate_compression_ratio_list_linear(input);
    generate_fix_compression_ratio_error_list(input);
}

void compress_fix_eb_best_cr(pressio_data *input, const double eb, const char *output_file) {

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

void compress_fix_psnr_best_cr(pressio_data *input, const double psnr, const char *output_file) {


}

void compress_fix_ssim_best_cr(pressio_data *input, const double ssim, const char *output_file) {

}

void compress_fix_cr_best_eb(pressio_data *input, const double cr, const char *output_file) {

    double best_eb = -1;
    std::string best_compressor = "null";

    for (const auto &_ : compressor_list) {

        double l = -1, r = 5, mid;
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

    compress_to_destination(best_compressor, input, best_eb, output_file);
}

void compress_fix_cr_best_psnr(pressio_data *input, const double cr, const char *output_file) {

}

void compress_fix_cr_best_ssim(pressio_data *input, const double cr, const char *output_file) {

}


double get_surrogate_cr(std::string _, pressio_data *input, const double eb) {

    if (is_reversed_dimension[_]) {
        pressio_data_reshape(input, reversed_dims.size(), reversed_dims.data());
    }
    else {
        pressio_data_reshape(input, dims.size(), dims.data());
    }

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
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double cr;
    assert_defined((_ + "_surrogate:cr").c_str(), &cr);
    printf("%s : Compression ratio: %lf with error bound %lf\n", compressor_string.c_str(), cr, eb);
    return cr;
}

double get_surrogate_psnr(const std::string _, pressio_data *input, const double eb) {

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
            throw std::runtime_error(key + " was not set"s);
        }
    };

    double psnr;
    assert_defined((_ + "_surrogate:psnr").c_str(), &psnr);
    // printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
    return psnr;
}

double get_real_cr(const std::string _, pressio_data *input, const double eb) {

    pressio_data_reshape(input, dims.size(), dims.data());

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

    double cr = 1. * input->size_in_bytes() / compressed.size_in_bytes();
    return cr;


    // auto metrics = compressor->get_metrics_results();
    // auto assert_defined = [&](const char* key, auto value){
    //     if(metrics.get(key, value)!= pressio_options_key_set) {
    //         throw std::runtime_error(key + " was not set"s);
    //     }
    // };

    // double cr;
    // // assert_defined("error_stat:psnr", &cr);
    // assert_defined("size:compression_ratio", &cr);
    // return cr;
}

void compress_to_destination(const std::string _, pressio_data *input, const double eb, const char *output_file) {

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
            throw std::runtime_error(key + " was not set"s);
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
//             throw std::runtime_error(key + " was not set"s);
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
//             throw std::runtime_error(key + " was not set"s);
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
//             throw std::runtime_error(key + " was not set"s);
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
//             throw std::runtime_error(key + " was not set"s);
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
//             throw std::runtime_error(key + " was not set"s);
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
//             throw std::runtime_error(key + " was not set"s);
//         }
//     };

//     double cr;
//     assert_defined("size:compression_ratio", &cr);
//     printf("%s : Compression ratio: %lf\n", compressor_string.c_str(), cr);
// }