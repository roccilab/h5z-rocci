#pragma once

#include "SZ3/utils/Config.hpp"
#include <stdexcept>

namespace SZ3_CR{double estimate_cr_float(SZ3::Config conf, float *data, double abs, int stride, uint N);};

namespace SZ3_PSNR{double estimate_psnr_float(SZ3::Config conf, float *data, double abs, int stride, uint N);};

namespace SZ3_SSIM{double estimate_ssim_float(SZ3::Config conf, float *data, double abs, int stride, int blocksize, uint N);};