#include <libpressio_ext/cpp/libpressio.h>

double zfp_estimate_cr_float(const struct pressio_data* input, double eb, double sample_rate);
double zfp_estimate_psnr_float(const struct pressio_data* input, double eb, double sample_rate);
double zfp_estimate_ssim_float(const struct pressio_data* input, double eb, double sample_rate);