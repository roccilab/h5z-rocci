#include <sstream>
#include <std_compat/memory.h>
#include <libpressio_ext/cpp/libpressio.h>
#include <SZ3/estimate_metric.h>


namespace libpressio { namespace sz3_surrogate {

    #define CR_METRIC 0
    #define PSNR_METRIC 1
    #define SSIM_METRIC 2


  const std::map<std::string,int> METRIC_MODES {
    {"cr", CR_METRIC},
    {"psnr", PSNR_METRIC},
    {"ssim", SSIM_METRIC},
  };

    struct iless {
        bool operator()(std::string lhs, std::string rhs) const {
            std::transform(std::begin(lhs), std::end(lhs), std::begin(lhs), [](unsigned char c){return std::tolower(c);});
            std::transform(std::begin(rhs), std::end(rhs), std::begin(rhs), [](unsigned char c){return std::tolower(c);});
            return lhs < rhs;
        }
    };  
  
struct sz3_option_maps {
  std::map<std::string, SZ3::EB, iless> error_bounds;
  std::map<std::string, SZ3::ALGO, iless> algo;
  std::map<std::string, SZ3::INTERP_ALGO, iless> interp_algo;
  sz3_option_maps() {
    for (size_t i = 0; i < std::size(SZ3::EB_STR); ++i) {
      error_bounds[SZ3::EB_STR[i]] = SZ3::EB_OPTIONS[i];
    }
    for (size_t i = 0; i < std::size(SZ3::ALGO_STR); ++i) {
      algo[SZ3::ALGO_STR[i]] = SZ3::ALGO_OPTIONS[i];
    }
    for (size_t i = 0; i < std::size(SZ3::INTERP_ALGO_STR); ++i) {
      interp_algo[SZ3::INTERP_ALGO_STR[i]] = SZ3::INTERP_ALGO_OPTIONS[i];
    }
  }
};
sz3_option_maps const& sz3_options() {
  static sz3_option_maps maps;
  return maps;
}
template <class K, class V, class Sort>
std::vector<K> keys(std::map<K,V,Sort> const& map) {
  std::vector<K> k;
  k.reserve(map.size());
  for (auto const& i : map) {
    k.emplace_back(i.first);
  }
  return k;
}

class sz3_surrogate_plugin : public libpressio_compressor_plugin {
public:
  sz3_surrogate_plugin() {
    config.absErrorBound = 1e-6;
    config.errorBoundMode = SZ3::EB_ABS;
  }

  struct pressio_options get_options_impl() const override
  {
    struct pressio_options options;
    if(config.errorBoundMode == SZ3::EB_ABS) {
      set(options, "pressio:abs", config.absErrorBound);
    } else {
      set_type(options, "pressio:abs", pressio_option_double_type);
    }
    if(config.errorBoundMode == SZ3::EB_REL) {
      set(options, "pressio:rel", config.relErrorBound);
    } else {
      set_type(options, "pressio:rel", pressio_option_double_type);
    }
    set(options, "pressio:nthreads", nthreads);
    set(options, "sz3_surrogate:abs_error_bound", config.absErrorBound);
    set(options, "sz3_surrogate:rel_error_bound", config.relErrorBound);
    set(options, "sz3_surrogate:error_bound_mode", config.errorBoundMode);
    set(options, "sz3_surrogate:quant_bin_size", config.quantbinCnt);
    set(options, "sz3_surrogate:blockSize", blockSize);
    set(options, "sz3_surrogate:samplingStride", sampleStride);
    set_type(options, "sz3_surrogate:metric", pressio_option_charptr_type);
    set_type(options, "sz3_surrogate:error_bound_mode_str", pressio_option_charptr_type);
    return options;
  }

  struct pressio_options get_configuration_impl() const override
  {
    struct pressio_options options;
    set(options, "pressio:thread_safe", pressio_thread_safety_multiple);
    set(options, "pressio:stability", "experimental");
    set(options, "sz3_surrogate:error_bound_mode_str", keys(sz3_options().error_bounds));
    set(options, "sz3_surrogate:metric", keys(METRIC_MODES));

    set(options, "predictors:error_dependent", std::vector<std::string>{
        "pressio:abs",
        "sz3_surrogate:abs_error_bound",
        "sz3_surrogate:rel_error_bound",
        "sz3_surrogate:err_bound_mode",
        "sz3_surrogate:err_bound_mode_str",
        "sz3_surrogate:blockSize",
        "sz3_surrogate:samplingStride",
        "sz3_surrogate:metric",
    });
    set(options, "predictors:error_agnostic", std::vector<std::string>{
        "pressio:abs",
        "sz3_surrogate:abs_error_bound",
        "sz3_surrogate:rel_error_bound",
        "sz3_surrogate:err_bound_mode",
        "sz3_surrogate:err_bound_mode_str",
        "sz3_surrogate:blockSize",
        "sz3_surrogate:samplingStride",
        "sz3_surrogate:metric",
    });

    return options;
  }

  struct pressio_options get_documentation_impl() const override
  {
    struct pressio_options options;
    set(options, "pressio:description", R"(A surrogate for SZ3)");
    set(options, "sz3_surrogate:abs_error_bound", "absolute error bound");
    set(options, "sz3_surrogate:rel_error_bound", "value range relative error bound");
    set(options, "sz3_surrogate:error_bound_mode", "error bound mode to apply");
    set(options, "sz3_surrogate:quant_bin_size", "number of quantization bins");
    set(options, "sz3_surrogate:error_bound_mode_str", "error bound");
    set(options, "sz3_surrogate:blockSize", "sample block size (only for ssim estimation)");
    set(options, "sz3_surrogate:samplingStride", "sample stride in interpolation estimator");
    set(options, "sz3_surrogate:metric", "metric to estimate, one of {cr, psnr, ssim}");
    return options;
  }


  int set_options_impl(struct pressio_options const& options) override
  {
    if(get(options, "pressio:abs", &config.absErrorBound) == pressio_options_key_set) {
      config.errorBoundMode = SZ3::EB_ABS;
    } 
    if(get(options, "pressio:rel", &config.relErrorBound) == pressio_options_key_set) {
      config.errorBoundMode = SZ3::EB_REL;
    } 
    uint32_t tmp_nthreads;
    if(get(options, "pressio:nthreads", &tmp_nthreads) == pressio_options_key_set) {
        if(tmp_nthreads > 1) {
            nthreads = tmp_nthreads;
            config.openmp = true;
        } else if(tmp_nthreads == 1) {
            nthreads = tmp_nthreads;
            config.openmp = false;
        } else {
            return set_error(1, "unsupported nthreads");
        }
    }
    get(options, "sz3_surrogate:abs_error_bound", &config.absErrorBound);
    get(options, "sz3_surrogate:rel_error_bound", &config.relErrorBound);
    get(options, "sz3_surrogate:error_bound_mode", &config.errorBoundMode);
    get(options, "sz3_surrogate:quant_bin_size", &config.quantbinCnt);
    get(options, "sz3_surrogate:blockSize", &blockSize);
    get(options, "sz3_surrogate:samplingStride", &sampleStride);
    std::string tmp;
    try {
      if(get(options, "sz3_surrogate:error_bound_mode_str", &tmp) == pressio_options_key_set) {
        config.errorBoundMode = sz3_options().error_bounds.at(tmp);
      }
    } catch(std::out_of_range const& ex) {
      return set_error(1, ex.what());
    }
    try {
      if(get(options, "sz3_surrogate:metric", &tmp) == pressio_options_key_set) {
        metric_setting = METRIC_MODES.at(tmp);
      }
    } catch(std::out_of_range const& ex) {
      return set_error(1, ex.what());
    }
    return 0;
  }

  int compress_impl(const pressio_data* input,
                    struct pressio_data* output) override
  {
    auto dims = input->dimensions();
    config.setDims(dims.begin(), dims.end());
    
    // avoid zero eb
    double eb = config.absErrorBound != 0 ? config.absErrorBound : 1e-16;
    if(config.errorBoundMode == SZ3::EB_REL) {
        printf("sz3_surrogate does not support a relative error bound\n");
        exit(1);
    }
    float* data = (float*) input->data();
    if (metric_setting == CR_METRIC) {
        estCR = estimate_cr_float(config, data, eb, sampleStride, input->num_dimensions());
    } else if (metric_setting == PSNR_METRIC) {
        estPSNR = estimate_psnr_float(config, data, eb, sampleStride, input->num_dimensions());
    } else if (metric_setting == SSIM_METRIC) {
        estSSIM = estimate_ssim_float(config, data, eb, sampleStride, blockSize, input->num_dimensions());
    } else {
        return set_error(1, "invalid metric string, expected one of {cr, psnr, ssim}");
    }

    return 0;

  }

  int decompress_impl(const pressio_data* real_input,
                      struct pressio_data* output) override
  {
        return 0;
  }

  int major_version() const override { return 0; }
  int minor_version() const override { return 0; }
  int patch_version() const override { return 0; }
  int revision_version() const { return 0; }
  const char* version() const override { 
    static std::string version_str = [this]{
      std::stringstream ss;
      ss << major_version() << '.' << minor_version() << '.' << patch_version() << '.' << revision_version();
      return ss.str();
    }();
    return version_str.c_str();
  }
  const char* prefix() const override { return "sz3_surrogate"; }

 pressio_options get_metrics_results_impl() const override {
    struct pressio_options options;
    set(options, "sz3_surrogate:cr", estCR);
    set(options, "sz3_surrogate:psnr", estPSNR);
    set(options, "sz3_surrogate:ssim", estSSIM);

    return options;
}

  std::shared_ptr<libpressio_compressor_plugin> clone() override
  {
    return compat::make_unique<sz3_surrogate_plugin>(*this);
  }

  uint32_t nthreads = 1;
  SZ3::Config config{};
  int sampleStride = 5;
  int blockSize = 32;
  int32_t metric_setting = CR_METRIC;
  double estCR = -1;
  double estPSNR = -1;
  double estSSIM = -1;
};

static pressio_register compressor_many_fields_plugin(compressor_plugins(), "sz3_surrogate", []() {
  return compat::make_unique<sz3_surrogate_plugin>();
});

} }