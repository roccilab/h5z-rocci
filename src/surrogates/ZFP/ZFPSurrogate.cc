#include <vector>
#include <memory>
#include <ZFP/estimate_metric.h>
#include <map>
#include <sstream>
#include <std_compat/memory.h>
#include <libpressio_ext/cpp/libpressio.h>
#include "pressio_options.h"
#include "pressio_data.h"
#include "pressio_compressor.h"
#include "pressio_version.h"
#include "std_compat/memory.h"
#include "std_compat/utility.h"

namespace libpressio { namespace zfp_surrogate_plugin {
  
  #define CR_METRIC 0
  #define PSNR_METRIC 1
  #define SSIM_METRIC 2

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

class zfp_surrogate_plugin: public libpressio_compressor_plugin {
  public:

    struct pressio_options get_options_impl() const override {
      struct pressio_options options;
      set_type(options, "pressio:abs", pressio_option_double_type);
      set_type(options, "zfp_surrogate:accuracy", pressio_option_double_type);
      set_type(options, "zfp_surrogate:sample_rate", pressio_option_double_type);
      set_type(options, "zfp_surrogate:metric", pressio_option_charptr_type);
      return options;
    }

    struct pressio_options get_documentation_impl() const override {
      struct pressio_options options;
      set(options, "pressio:description", R"(Surrogate for ZFP)");
      set(options, "zfp_surrogate:accuracy", "absolute error tolerance for fixed-accuracy mode");
      set(options, "zfp_surrogate:sample_rate", "Sample rate in [0,1] for surrogate");
      set(options, "zfp_surrogate:metric", "Metric to estimate, one of {cr, psnr, ssim}");
      return options;
    }

    struct pressio_options get_configuration_impl() const override {
      struct pressio_options options;
      set(options, "pressio:thread_safe", pressio_thread_safety_multiple);

      set(options, "predictors:error_dependent", std::vector<std::string>{
              "pressio:abs",
              "zfp_surrogate:accuracy",
              "zfp_surrogate:sample_rate",
              });
      set(options, "predictors:error_agnostic", std::vector<std::string>{
              "pressio:abs",
              "zfp_surrogate:accuracy",
              "zfp_surrogate:sample_rate",
              });
    return options;
    }

    int set_options_impl(struct pressio_options const& options) override {
      double tolerance; 
      if(get(options, "pressio:abs", &tolerance) == pressio_options_key_set) {
        accuracy = tolerance;
      } else if (get(options, "zfp_surrogate:accuracy", &tolerance) == pressio_options_key_set) {
        accuracy = tolerance;
      } 

      get(options, "zfp_surrogate:sample_rate", &sample_rate);

      std::string tmp;
      if(get(options, "zfp_surrogate:metric", &tmp) == pressio_options_key_set) {
        metric_setting = METRIC_MODES.at(tmp);
      }
      return 0;
    }

    int compress_impl(const pressio_data *input, struct pressio_data* output) override {

        if (metric_setting == CR_METRIC) {
            estCR = zfp_estimate_cr_float(input, accuracy, sample_rate);
        } else if (metric_setting == PSNR_METRIC) {
            estPSNR = zfp_estimate_psnr_float(input, accuracy, sample_rate);
        } else if (metric_setting == SSIM_METRIC) {
            estSSIM = zfp_estimate_ssim_float(input, accuracy, sample_rate);
        } else {
            return set_error(1, "invalid metric string, expected one of {cr, psnr, ssim}");
        }
        
        return 0;
    }

    int decompress_impl(const pressio_data *input, struct pressio_data* output) override {
      
      return 0;
    }

    pressio_options get_metrics_results_impl() const override {
        struct pressio_options options;
        set(options, "zfp_surrogate:cr", estCR);
        set(options, "zfp_surrogate:psnr", estPSNR);
        set(options, "zfp_surrogate:ssim", estSSIM);

        return options;
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

    const char* prefix() const override {
      return "zfp_surrogate";
    }

    std::shared_ptr<libpressio_compressor_plugin> clone() override{
      return compat::make_unique<zfp_surrogate_plugin>(*this);
    }

    double accuracy;
    double sample_rate = 0.2;
    double estCR = -1, estPSNR = -1, estSSIM = -1;
    int32_t metric_setting;
};

static pressio_register compressor_many_fields_plugin(compressor_plugins(), "zfp_surrogate", [](){ return compat::make_unique<zfp_surrogate_plugin>(); });

} }