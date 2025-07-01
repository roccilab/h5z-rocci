#pragma once

static const int cmp_tblock_size_f32 = 32;
static const int dec_tblock_size_f32 = 32;
static const int cmp_chunk_f32 = 8192;
static const int dec_chunk_f32 = 8192;

double SZp_estimate_compress_hostptr_f32(float* oriData, size_t nbEle, float errorBound, size_t sample_stride);
double SZp_estimate_psnr_hostptr_f32(float* oriData, size_t nbEle, float errorBound, size_t sample_stride);
double SZp_estimate_ssim_hostptr_f32(float* oriData, size_t dims[3], size_t n_dim, size_t nbEle, float errorBound, size_t sample_stride);

