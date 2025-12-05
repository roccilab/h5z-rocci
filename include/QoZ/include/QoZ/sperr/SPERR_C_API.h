/*
 * This header provides C API for SPERR.
 * This API is supposed to be used in non-C++ (e.g., C, Fortran, Python) projects.
 */

#ifndef SPERR_C_API_H
#define SPERR_C_API_H

#include "SPERR2D_Compressor.h"
#include "SPERR2D_Decompressor.h"

#include "SPERR3D_OMP_C.h"
#include "SPERR3D_OMP_D.h"
#include <stddef.h> /* for size_t */
#include <stdint.h> /* for fixed-width integers */

#ifdef __cplusplus
namespace C_API {
extern "C" {
#endif

/*
 * The memory management is a little tricy. The following requirement applies to the output
 * buffer `dst` of all the C API functions.
 *
 * The output is stored in `dst`, which is a pointer pointing to another pointer
 * held by the caller.  The other pointer should be NULL; otherwise this function will fail!
 * Upon success, `dst` will contain a buffer of length `dst_len` in case of compression,
 * and (dimx x dimy x dimz) of floats (or doubles) in case of decompression.
 * The caller of this function is responsible of free'ing `dst` using free().
 *
 */

/*
 * Compress a a 2D slice targetting different quality controls (modes):
 *   mode == 1 --> fixed bit-per-pixel (BPP)
 *   mode == 2 --> fixed peak signal-to-noise ratio (PSNR)
 *   mode == 3 --> fixed point-wise error (PWE)
 *
 * Return value meanings:
 * 0: success
 * 1: `dst` is not pointing to a NULL pointer!
 * 2: `mode` or `quality` isn't valid
 * 3: `is_float` value not supported
 *-1: other error
 */
int sperr_comp_2d(
    const void* src,  /* Input: buffer that contains a 2D slice */
    int32_t is_float, /* Input: input buffer type: 1 == float, 0 == double */
    size_t dimx,      /* Input: X (fastest-varying) dimension */
    size_t dimy,      /* Input: Y (slowest-varying) dimension */
    int32_t mode,     /* Input: compression mode to use */
    double quality,   /* Input: target quality */
    void** dst,       /* Output: buffer for the output bitstream, allocated by this function */
    size_t* dst_len); /* Output: length of `dst` in byte */

/*
 * Compress a a 3D volume targetting different quality controls (modes):
 *   mode == 1 --> fixed bit-per-pixel (BPP)
 *   mode == 2 --> fixed peak signal-to-noise ratio (PSNR)
 *   mode == 3 --> fixed point-wise error (PWE)
 *
 * Return value meanings:
 * 0: success
 * 1: `dst` is not pointing to a NULL pointer!
 * 2: `mode` or `quality` isn't valid
 * 3: `is_float` value not supported
 *-1: other error
 */
int sperr_comp_3d(
    const void* src,  /* Input: buffer that contains a 3D volume */
    int32_t is_float, /* Input: input buffer type: 1 == float, 0 = double */
    size_t dimx,      /* Input: X (fastest-varying) dimension */
    size_t dimy,      /* Input: Y dimension */
    size_t dimz,      /* Input: Z (slowest-varying) dimension */
    size_t chunk_x,   /* Input: preferred chunk dimension in X */
    size_t chunk_y,   /* Input: preferred chunk dimension in Y */
    size_t chunk_z,   /* Input: preferred chunk dimension in Z */
    int32_t mode,     /* Input: compression mode to use */
    double quality,   /* Input: target quality */
    size_t nthreads,  /* Input: number of OMP threads to use. 0 means using all threads. */
    void** dst,       /* Output: buffer for the output bitstream, allocated by this function */
    size_t* dst_len); /* Output: length of `dst` in byte */

/*
 * Decompress a 2D SPERR-compressed buffer.
 *
 * Return value meanings:
 * 0: success
 * 1: `dst` not pointing to a NULL pointer!
 * 2: `output_float` value not supported
 * -1: other error
 */
int sperr_decomp_2d(
    const void* src,      /* Input: buffer that contains a compressed bitstream */
    size_t src_len,       /* Input: length of the input bitstream in byte */
    int32_t output_float, /* Input: output data type: 1 == float, 0 == double */
    size_t* dimx,         /* Output: X (fast-varying) dimension */
    size_t* dimy,         /* Output: Y (slowest-varying) dimension */
    void** dst);          /* Output: buffer for the output 2D slice, allocated by this function */

/*
 * Decompress a 3D SPERR-compressed buffer.
 *
 * Return value meanings:
 * 0: success
 * 1: `dst` not pointing to a NULL pointer!
 * 2: `output_float` value not supported
 * -1: other error
 */
int sperr_decomp_3d(
    const void* src,      /* Input: buffer that contains a compressed bitstream */
    size_t src_len,       /* Input: length of the input bitstream in byte */
    int32_t output_float, /* Input: output data type: 1 == float, 0 == double */
    size_t nthreads,      /* Input: number of OMP threads to use. 0 means using all threads. */
    size_t* dimx,         /* Output: X (fast-varying) dimension */
    size_t* dimy,         /* Output: Y dimension */
    size_t* dimz,         /* Output: Z (slowest-varying) dimension */
    void** dst);          /* Output: buffer for the output 3D slice, allocated by this function */

/*
 * Decompress a 2D or 3D SPERR bitstream to a USER-ALLOCATED memory buffer. If you don't know
 * how much memory to allocate, use functions `sperr_decomp_2d()` or `sperr_decomp_3d()` instead.
 * Note 1: If not enough memory was allocated at `dst`, segment faults may occur.
 * Note 2: `nthreads` won't have any effect in 2D decompression.
 *
 * Return value meanings:
 * 0: success
 * 1: `output_float` value not supported
 *-1: other error
 */
int sperr_decomp_user_mem(
    const void* src,      /* Input: buffer that contains a compressed bitstream */
    size_t src_len,       /* Input: length of the input bitstream in byte */
    int32_t output_float, /* Input: output data type: 1 == float, 0 == double */
    size_t nthreads,      /* Input: number of OMP threads to use. 0 means using all threads. */
    void* dst);           /* Output: decompressed slice or volume in memory allocated by caller */

/*
 * Given a SPERR bitstream, parse the header and retrieve various types of information.
 *
 * Return value meanings:
 * 0: success
 * 1: parsing error occured
 */
void sperr_parse_header(
    const void* ptr,        /* Input: the bitstream to parse */
    int32_t* version_major, /* Output: major version number */
    int32_t* zstd_applied,  /* Output: if ZSTD applied (0 == no ZSTD; 1 == ZSTD applied) */
    int32_t* is_3d,         /* Output: 3D volume or 2D slice (0 == 2D, 1 == 3D) */
    int32_t* orig_is_float, /* Output: Precision of the original input (0 == double, 1 == float) */
    uint32_t* dim_x,        /* Output: X dimension */
    uint32_t* dim_y,        /* Output: Y dimension */
    uint32_t* dim_z);       /* Output: Z dimension (undefined value for 2D slices) */

#ifdef __cplusplus
} /* end of extern "C" */
} /* end of namespace C_API */
#endif



int C_API::sperr_comp_2d(const void* src,
                         int32_t is_float,
                         size_t dimx,
                         size_t dimy,
                         int32_t mode,
                         double quality,
                         void** dst,
                         size_t* dst_len)
{
  // Examine if `dst` is pointing to a NULL pointer
  if (*dst != NULL)
    return 1;

  // Examine if `mode` and `quality` are valid
  if (mode < 1 || mode > 3 || quality <= 0.0)
    return 2;

  // Setup the compressor
  const auto total_vals = dimx * dimy;
  auto compressor = SPERR2D_Compressor();
  auto rtn = sperr::RTNType::Good;
  switch (is_float) {
    case 0:  // double
      rtn = compressor.copy_data(static_cast<const double*>(src), total_vals, {dimx, dimy, 1});
      break;
    case 1:  // float
      rtn = compressor.copy_data(static_cast<const float*>(src), total_vals, {dimx, dimy, 1});
      break;
    default:
      rtn = RTNType::Error;
  }
  if (rtn != RTNType::Good)
    return -1;

  // Specify a particular compression mode.
  switch (mode) {
    case 1:
      compressor.set_target_bpp(quality);
      break;
    case 2:
      compressor.set_target_psnr(quality);
      break;
    case 3:
      compressor.set_target_pwe(quality);
      break;
    default:
      return 2;
  }

  // Do the actual compression work!
  rtn = compressor.compress();
  if (rtn != RTNType::Good)
    return -1;

  // Output the compressed bitstream
  const auto& stream = compressor.view_encoded_bitstream();
  if (stream.empty())
    return -1;
  *dst_len = stream.size();
  uint8_t* buf = (uint8_t*)std::malloc(stream.size());
  std::copy(stream.cbegin(), stream.cend(), buf);
  *dst = buf;

  return 0;
}

int C_API::sperr_comp_3d(const void* src,
                         int32_t is_float,
                         size_t dimx,
                         size_t dimy,
                         size_t dimz,
                         size_t chunk_x,
                         size_t chunk_y,
                         size_t chunk_z,
                         int32_t mode,
                         double quality,
                         size_t nthreads,
                         void** dst,
                         size_t* dst_len)
{
  // Examine if `dst` is pointing to a NULL pointer
  if (*dst != NULL)
    return 1;

  // Examine if `mode` and `quality` are valid
  if (mode < 1 || mode > 3 || quality <= 0.0)
    return 2;

  const auto vol_dims = sperr::dims_type{dimx, dimy, dimz};
  const auto chunk_dims = sperr::dims_type{chunk_x, chunk_y, chunk_z};

  // Setup the compressor
  const auto total_vals = dimx * dimy * dimz;
  auto compressor = SPERR3D_OMP_C();
  compressor.set_num_threads(nthreads);
  auto rtn = sperr::RTNType::Good;
  switch (is_float) {
    case 0:  // double
      rtn = compressor.copy_data(static_cast<const double*>(src), total_vals, vol_dims, chunk_dims);
      break;
    case 1:  // float
      rtn = compressor.copy_data(static_cast<const float*>(src), total_vals, vol_dims, chunk_dims);
      break;
    default:
      rtn = RTNType::Error;
  }
  if (rtn != RTNType::Good)
    return -1;

  // Specify a particular compression mode.
  switch (mode) {
    case 1:
      compressor.set_target_bpp(quality);
      break;
    case 2:
      compressor.set_target_psnr(quality);
      break;
    case 3:
      compressor.set_target_pwe(quality);
      break;
    default:
      return 2;
  }

  // Do the actual compression work
  rtn = compressor.compress();
  if (rtn != RTNType::Good)
    return -1;

  // Output the compressed bitstream
  const auto stream = compressor.get_encoded_bitstream();
  if (stream.empty())
    return -1;
  *dst_len = stream.size();
  uint8_t* buf = (uint8_t*)std::malloc(stream.size());
  std::copy(stream.cbegin(), stream.cend(), buf);
  *dst = buf;

  return 0;
}

int C_API::sperr_decomp_2d(const void* src,
                           size_t src_len,
                           int32_t output_float,
                           size_t* dimx,
                           size_t* dimy,
                           void** dst)
{
  // Examine if `dst` is pointing to a NULL pointer
  if (*dst != NULL)
    return 1;

  // Use a decompressor to decompress this bitstream
  auto decompressor = SPERR2D_Decompressor();
  auto rtn = decompressor.use_bitstream(src, src_len);
  if (rtn != RTNType::Good)
    return -1;
  rtn = decompressor.decompress();
  if (rtn != RTNType::Good)
    return -1;

  // Double check that the slice dimension is correct
  const auto& slice = decompressor.view_data();
  const auto dims = decompressor.get_dims();
  if (slice.size() != dims[0] * dims[1])
    return -1;
  *dimx = dims[0];
  *dimy = dims[1];

  // write out the 2D slice in double or float format
  switch (output_float) {
    case 0: {  // double
      double* buf = (double*)std::malloc(slice.size() * sizeof(double));
      std::copy(slice.cbegin(), slice.cend(), buf);
      *dst = buf;
      break;
    }
    case 1: {  // float
      float* buf = (float*)std::malloc(slice.size() * sizeof(float));
      std::copy(slice.cbegin(), slice.cend(), buf);
      *dst = buf;
      break;
    }
    default:
      return 2;
  }

  return 0;
}

int C_API::sperr_decomp_3d(const void* src,
                           size_t src_len,
                           int32_t output_float,
                           size_t nthreads,
                           size_t* dimx,
                           size_t* dimy,
                           size_t* dimz,
                           void** dst)
{
  // Examine if `dst` is pointing to a NULL pointer
  if (*dst != NULL)
    return 1;

  // Use a decompressor to decompress this bitstream
  auto decompressor = SPERR3D_OMP_D();
  decompressor.set_num_threads(nthreads);
  auto rtn = decompressor.use_bitstream(src, src_len);
  if (rtn != RTNType::Good)
    return -1;
  rtn = decompressor.decompress(src);
  if (rtn != RTNType::Good)
    return -1;

  // Double check that the volume dimension is correct
  const auto& vol = decompressor.view_data();
  const auto dims = decompressor.get_dims();
  if (vol.size() != dims[0] * dims[1] * dims[2])
    return -1;
  *dimx = dims[0];
  *dimy = dims[1];
  *dimz = dims[2];

  // write out the volume in double or float format
  switch (output_float) {
    case 0: {  // double
      double* buf = (double*)std::malloc(vol.size() * sizeof(double));
      std::copy(vol.cbegin(), vol.cend(), buf);
      *dst = buf;
      break;
    }
    case 1: {  // float
      float* buf = (float*)std::malloc(vol.size() * sizeof(float));
      std::copy(vol.cbegin(), vol.cend(), buf);
      *dst = buf;
      break;
    }
    default:
      return 2;
  }

  return 0;
}

int C_API::sperr_decomp_user_mem(const void* src,
                                 size_t src_len,
                                 int32_t output_float,
                                 size_t nthreads,
                                 void* dst)
{
  auto header = sperr::parse_header(src);
  if (output_float < 0 || output_float > 1)
    return 1;

  if (header.is_3d) {  // Decompress a 3D volume
    auto decompressor = SPERR3D_OMP_D();
    decompressor.set_num_threads(nthreads);
    auto rtn = decompressor.use_bitstream(src, src_len);
    if (rtn != RTNType::Good)
      return -1;
    rtn = decompressor.decompress(src);
    if (rtn != RTNType::Good)
      return -1;

    const auto& vol = decompressor.view_data();
    if (output_float == 0) {  // Output in double precision
      double* const ptr = static_cast<double*>(dst);
      std::copy(vol.cbegin(), vol.cend(), ptr);
    }
    else {  // Output in single precisioin
      float* const ptr = static_cast<float*>(dst);
      std::copy(vol.cbegin(), vol.cend(), ptr);
    }
  }
  else {  // Decompress a 2D slice
    auto decompressor = SPERR2D_Decompressor();
    auto rtn = decompressor.use_bitstream(src, src_len);
    if (rtn != RTNType::Good)
      return -1;
    rtn = decompressor.decompress();
    if (rtn != RTNType::Good)
      return -1;

    const auto& slice = decompressor.view_data();
    if (output_float == 0) {  // Output in double precision
      double* const ptr = static_cast<double*>(dst);
      std::copy(slice.cbegin(), slice.cend(), ptr);
    }
    else {  // Output in single precisioin
      float* const ptr = static_cast<float*>(dst);
      std::copy(slice.cbegin(), slice.cend(), ptr);
    }
  }

  return 0;
}

void C_API::sperr_parse_header(const void* ptr,
                               int32_t* version_major,
                               int32_t* zstd_applied,
                               int32_t* is_3d,
                               int32_t* orig_is_float,
                               uint32_t* dim_x,
                               uint32_t* dim_y,
                               uint32_t* dim_z)
{
  auto header = sperr::parse_header(ptr);

  *version_major = header.version_major;
  *zstd_applied = header.zstd_applied;
  *is_3d = header.is_3d;
  *orig_is_float = header.orig_is_float;
  *dim_x = header.vol_dims[0];
  *dim_y = header.vol_dims[1];
  *dim_z = header.vol_dims[2];
}



#endif
