#pragma once

#include "SZ3/utils/Config.hpp"
#include <stdexcept>

double estimate_cr_float(SZ3::Config conf, float *data, double abs, int stride, uint N);

double estimate_psnr_float(SZ3::Config conf, float *data, double abs, int stride, uint N);

double estimate_ssim_float(SZ3::Config conf, float *data, double abs, int stride, int blocksize, uint N);