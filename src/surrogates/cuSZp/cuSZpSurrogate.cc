#include <map>
#include <sstream>
#include <algorithm>
#include <std_compat/memory.h>
#include <libpressio_ext/cpp/libpressio.h>
#include <cuSZp/cuSZp_f32.h>


namespace libpressio { namespace cuszp_surrogate {

  #define ABS 0
  #define REL 1

  #define CR_METRIC 0
  #define PSNR_METRIC 1
  #define SSIM_METRIC 2

  const std::map<std::string,int> ERR_MODES {
    {"abs", ABS},
    {"rel", REL},
  };

  const std::map<std::string,int> METRIC_MODES {
    {"cr", CR_METRIC},
    {"psnr", PSNR_METRIC},
    {"ssim", SSIM_METRIC},
  };

  std::vector<std::string> keys(std::map<std::string, int> const& entries) {
    std::vector<std::string> k;
    k.reserve(entries.size());
    for (auto const& i : entries) {
      k.emplace_back(i.first);
    }
    return k;
  }

class cuszp_surrogate_plugin : public libpressio_compressor_plugin {
public:
  struct pressio_options get_options_impl() const override
  {
    struct pressio_options options;
    if(errBoundMode == ABS){
      set(options, "pressio:abs", absErrBound);
    } else {
      set_type(options, "pressio:abs", pressio_option_double_type);
    }
    if(errBoundMode == REL){
      set(options, "pressio:rel", relBoundRatio);
    } else {
      set_type(options, "pressio:rel", pressio_option_double_type);
    }

    set(options, "cuszp_surrogate:abs_err_bound", absErrBound);
    set(options, "cuszp_surrogate:rel_bound_ratio", relBoundRatio);
    set(options, "cuszp_surrogate:err_bound_mode", errBoundMode);
    set(options, "cuszp_surrogate:samplingStride", samplingStride);
    set_type(options, "cuszp_surrogate:metric", pressio_option_charptr_type);
    set_type(options, "cuszp_surrogate:err_bound_mode_str", pressio_option_charptr_type);
    return options;
  }

  struct pressio_options get_configuration_impl() const override
  {
    struct pressio_options options;
    set(options, "pressio:thread_safe", pressio_thread_safety_multiple);
    set(options, "pressio:stability", "experimental");
    set(options, "cuszp_surrogate:err_bound_mode_str", keys(ERR_MODES));
    set(options, "cuszp_surrogate:metric", keys(METRIC_MODES));
    


    set(options, "predictors:error_dependent", std::vector<std::string>{
    "pressio:abs",
    "cuszp_surrogate:abs_err_bound",
    "cuszp_surrogate:rel_bound_ratio",
    "cuszp_surrogate:err_bound_mode",
    "cuszp_surrogate:err_bound_mode_str",
    "cuszp_surrogate:samplingStride",
    "cuszp_surrogate:metric",
    });
      set(options, "predictors:error_agnostic", std::vector<std::string>{
    "pressio:abs",
    "cuszp_surrogate:abs_err_bound",
    "cuszp_surrogate:rel_bound_ratio",
    "cuszp_surrogate:err_bound_mode",
    "cuszp_surrogate:err_bound_mode_str",
    "cuszp_surrogate:samplingStride",
    "cuszp_surrogate:metric",
    });

    return options;
  }

  struct pressio_options get_documentation_impl() const override
  {
    struct pressio_options options;
    set(options, "pressio:description", R"(A surrogate for cuSZp)");
    set(options, "cuszp_surrogate:abs_err_bound", "absolute pointwise error bound");
    set(options, "cuszp_surrogate:rel_bound_ratio", "pointwise relative error bound error bound");
    set(options, "cuszp_surrogate:err_bound_mode", "error bound type");
    set(options, "cuszp_surrogate:err_bound_mode_str", "error bound type as a human readable string");
    set(options, "cuszp_surrogate:samplingStride", "sampling stride for estimation");
    set(options, "cuszp_surrogate:metric", "string describing metric to estimate - {'cr', 'psnr', 'ssim'}");
    
    return options;
  }


  int set_options_impl(struct pressio_options const& options) override
  {
    try {
      if(get(options, "pressio:abs", &absErrBound) == pressio_options_key_set) {
        errBoundMode = ABS;
      }
      if(get(options, "pressio:rel", &relBoundRatio) == pressio_options_key_set) {
        errBoundMode = REL;
      }

      get(options, "cuszp_surrogate:abs_err_bound", &absErrBound);
      get(options, "cuszp_surrogate:rel_bound_ratio", &relBoundRatio);
      get(options, "cuszp_surrogate:err_bound_mode", &errBoundMode);
      get(options, "cuszp_surrogate:samplingStride", &samplingStride);

      std::string tmp;
      if(get(options, "cuszp_surrogate:err_bound_mode_str", &tmp) == pressio_options_key_set) {
          errBoundMode = ERR_MODES.at(tmp);
      }
      if(get(options, "cuszp_surrogate:metric", &tmp) == pressio_options_key_set) {
        metric_setting = METRIC_MODES.at(tmp);
      }
    }
    catch (std::out_of_range const& ex) {
      return set_error(1, ex.what());
    }

    return 0;
  }

  int compress_impl(const pressio_data* input,
                    struct pressio_data* output) override
  {
    double eb = absErrBound;
    if(errBoundMode == REL) {
        eb = relBoundRatio;
        printf("cuszp_surrogate does not support a relative error bound\n");
        exit(1);
    }
    float* data = (float*) input->data();
    if (metric_setting == CR_METRIC) {
        estCR = SZp_estimate_compress_hostptr_f32(data, input->num_elements(), eb, samplingStride);
    } else if (metric_setting == PSNR_METRIC) {
        estPSNR = SZp_estimate_psnr_hostptr_f32(data, input->num_elements(), eb, samplingStride);
    } else if (metric_setting == SSIM_METRIC) {
        std::vector<size_t> input_dims(input->dimensions().rbegin(), input->dimensions().rend());
        estSSIM = SZp_estimate_ssim_hostptr_f32(data, input_dims.data(), input->num_dimensions(), input->num_elements(), eb, samplingStride);
    } else {
        return set_error(1, "invalid metric string, expected one of {cr, psnr, ssim}");
    }

    return 0;

  }

  int decompress_impl(const pressio_data* input,
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
  const char* prefix() const override { return "cuszp_surrogate"; }

  pressio_options get_metrics_results_impl() const override {
    struct pressio_options options;
    set(options, "cuszp_surrogate:cr", estCR);
    set(options, "cuszp_surrogate:psnr", estPSNR);
    set(options, "cuszp_surrogate:ssim", estSSIM);

    return options;
  }

  std::shared_ptr<libpressio_compressor_plugin> clone() override
  {
    return compat::make_unique<cuszp_surrogate_plugin>(*this);
  }

private:

  int32_t errBoundMode = ABS;
  double absErrBound = 1e-4;
  double relBoundRatio = 0;
  int32_t samplingStride = 5;
  int32_t metric_setting = CR_METRIC;
  double estCR = -1;
  double estPSNR = -1;
  double estSSIM = -1;
};

static pressio_register compressor_many_fields_plugin(compressor_plugins(), "cuszp_surrogate", []() {
  return compat::make_unique<cuszp_surrogate_plugin>();
});

} }