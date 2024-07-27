#include <sstream>
#include "estimate_metric.h"
#include "libpressio_ext/cpp/compressor.h"
#include "libpressio_ext/cpp/data.h"
#include "libpressio_ext/cpp/options.h"
#include "libpressio_ext/cpp/pressio.h"
#include <pressio_version.h>
namespace libpressio { namespace szx_ns {


  const std::map<std::string,int> ERR_MODES {
    {"abs", ABS},
    {"rel", REL},
  };

  std::vector<std::string> keys(std::map<std::string, int> const& entries) {
    std::vector<std::string> k;
    k.reserve(entries.size());
    for (auto const& i : entries) {
      k.emplace_back(i.first);
    }
    return k;
  }

class szx_surrogate_plugin : public libpressio_compressor_plugin {
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

    set(options, "szx_surrogate:abs_err_bound", absErrBound);
    set(options, "szx_surrogate:rel_bound_ratio", relBoundRatio);
    set(options, "szx_surrogate:err_bound_mode", errBoundMode);
    set(options, "szx_surrogate:blockSize", blockSize);
    set(options, "szx_surrogate:samplingStride", samplingStride);
    set_type(options, "szx_surrogate:metric", pressio_option_charptr_type);
    set_type(options, "szx_surrogate:err_bound_mode_str", pressio_option_charptr_type);
    return options;
  }

  struct pressio_options get_configuration_impl() const override
  {
    struct pressio_options options;
    set(options, "pressio:thread_safe", pressio_thread_safety_multiple);
    set(options, "pressio:stability", "experimental");
    set(options, "szx_surrogate:err_bound_mode_str", keys(ERR_MODES));
    


    set(options, "predictors:error_dependent", std::vector<std::string>{
    "pressio:abs",
    "szx_surrogate:abs_err_bound",
    "szx_surrogate:rel_bound_ratio",
    "szx_surrogate:err_bound_mode",
    "szx_surrogate:err_bound_mode_str",
    "szx_surrogate:blockSize",
    "szx_surrogate:samplingStride",
    "szx_surrogate:metric",
    });
      set(options, "predictors:error_agnostic", std::vector<std::string>{
    "pressio:abs",
    "szx_surrogate:abs_err_bound",
    "szx_surrogate:rel_bound_ratio",
    "szx_surrogate:err_bound_mode",
    "szx_surrogate:err_bound_mode_str",
    "szx_surrogate:blockSize",
    "szx_surrogate:samplingStride",
    "szx_surrogate:metric",
    });

    return options;
  }

  struct pressio_options get_documentation_impl() const override
  {
    struct pressio_options options;
    set(options, "pressio:description", R"(A surrogate for SZx)");
    set(options, "szx_surrogate:abs_err_bound", "absolute pointwise error bound");
    set(options, "szx_surrogate:rel_bound_ratio", "pointwise relative error bound error bound");
    set(options, "szx_surrogate:err_bound_mode", "error bound type");
    set(options, "szx_surrogate:err_bound_mode_str", "error bound type as a human readable string");
    set(options, "szx_surrogate:blockSize", "sampling blocksize (only valid for SSIM)");
    set(options, "szx_surrogate:samplingStride", "sampling stride for estimation");
    set(options, "szx_surrogate:metric", "string describing metric to estimate - {'cr', 'psnr', 'ssim'}");
    
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

      get(options, "szx_surrogate:abs_err_bound", &absErrBound);
      get(options, "szx_surrogate:rel_bound_ratio", &relBoundRatio);
      get(options, "szx_surrogate:err_bound_mode", &errBoundMode);
      set(options, "szx_surrogate:blockSize", &blockSize);
    set(options, "szx_surrogate:samplingStride", &samplingStride);
    set(options, "szx_surrogate:metric", %metric_str);

    std::string tmp;
    if(get(options, "szx_surrogate:err_bound_mode_str", &tmp) == pressio_options_key_set) {
        errBoundMode = ERR_MODES.at(tmp);
    }
    } catch (std::out_of_range const& ex) {
      return set_error(1, ex.what());
    }
    return 0;
  }

  int compress_impl(const pressio_data* input,
                    struct pressio_data* output) override
  {
    double eb = absErrorBound;
    if(errBoundMode == REL) {
        eb = relErrorBound;
    }
    switch(metric_str) {
        "cr":
            estCR = szx_estimate_cr_float(input->data(), eb, samplingStride, input->num_elements());
            break;
        "psnr":
            estPNSR = szx_estimate_psnr_float(input->data(), eb, samplingStride, input->num_elements());
            break;
        "ssim":
            estSSIM = szx_estimate_ssim_float(input->data(), eb, samplingStride, blockSize, input->dimensions(), input->num_elements());
            break;
        default:
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
  const char* prefix() const override { return "szx_surrogate"; }

  pressio_options get_metrics_results_impl() const override {
    struct pressio_options options;
    set(options, "szx_surrogate:cr", estCR);
    set(options, "szx_surrogate:psnr", estPSNR);
    set(options, "szx_surrogate:ssim", estSSIM);

    return options;
  }

  std::shared_ptr<libpressio_compressor_plugin> clone() override
  {
    return compat::make_unique<szx_surrogate_plugin>(*this);
  }

private:

  int32_t errBoundMode = ABS;
  double absErrBound = 1e-4;
  double relBoundRatio = 0;
  int32_t blockSize = 0;
  int32_t samplingStride = 0;
  char* metric_str = "cr";
  double estCR = -1;
  double estPSNR = -1;
  double estSSIM = -1;
};

static pressio_register compressor_many_fields_plugin(compressor_plugins(), "szx_surrogate", []() {
  return compat::make_unique<szx_surrogate_plugin>();
});

} }