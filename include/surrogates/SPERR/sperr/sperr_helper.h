#ifndef SPERR_HELPER_H
#define SPERR_HELPER_H

//
// We put common functions that are used across the project here.
//

#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // fixed width integers
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <utility>  // std::pair
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>

#ifdef USE_OMP
#include <omp.h>
#endif
//#ifndef USE_VANILLA_CONFIG
//#include "SperrConfig.h"
//#endif

namespace sperr {

using std::size_t;  // Seems most appropriate

// Shortcut for the maximum values
constexpr auto max_size = std::numeric_limits<size_t>::max();
constexpr auto max_d = std::numeric_limits<double>::max();

//
// A few shortcuts
//
using vecd_type = std::vector<double>;
using vec8_type = std::vector<uint8_t>;
using dims_type = std::array<size_t, 3>;

//
// Helper classes
//
enum class SigType : unsigned char { Insig, Sig, NewlySig, Dunno, Garbage };

enum class SetType : unsigned char { TypeS, TypeI, Garbage };

// Return Type
enum class RTNType {
  Good = 0,
  WrongDims,
  BitstreamWrongLen,
  IOError,
  InvalidParam,
  QzLevelTooBig,  // a very specific type of invalid param
  EmptyStream,    // a condition but not sure if it's an error
  BitBudgetMet,
  VersionMismatch,
  ZSTDMismatch,
  ZSTDError,
  SliceVolumeMismatch,
  QzModeMismatch,
  SetBPPBeforeDims,
  DataRangeNotSet,
  CompModeUnknown,
  CustomFilterMissing,
  CustomFilterError,
  Error
};

// Compression Mode
enum class CompMode { FixedSize, FixedQz, FixedPSNR, FixedPWE, Unknown };

//
// Helper functions
//
// Given a certain length, how many transforms to be performed?
auto num_of_xforms(size_t len) -> size_t;

// How many partition operation could we perform given a length?
// Length 0 and 1 can do 0 partitions; len=2 can do 1; len=3 can do 2, len=4 can
// do 2, etc.
auto num_of_partitions(size_t len) -> size_t;

// Determine the approximation and detail signal length at a certain
// transformation level lev: 0 <= lev < num_of_xforms.
// It puts the approximation and detail length as the 1st and 2nd
// element of the return array.
auto calc_approx_detail_len(size_t orig_len, size_t lev) -> std::array<size_t, 2>;

// Pack and unpack booleans to array of chars.
// When packing, the caller should make sure the number of booleans is a
// multiplier of 8. It optionally takes in an offset that specifies where to
// start writing/reading the char array.
//
// Note 1: unpack_booleans() takes a raw pointer because it accesses memory
// provided by others, and others most likely provide it by raw pointers.
// Note 2: these two methods only work on little endian machines.
// Note 3: the caller should have already allocated enough space for `dest`.
auto pack_booleans(std::vector<uint8_t>& dest, const std::vector<bool>& src, size_t dest_offset = 0)
    -> RTNType;
auto unpack_booleans(std::vector<bool>& dest,
                     const void* src,
                     size_t src_len,
                     size_t src_offset = 0) -> RTNType;

// Pack and unpack exactly 8 booleans to/from a single byte
// Note: memory for the 8 booleans should already be allocated!
// Note: these two methods only work on little endian machines.
auto pack_8_booleans(std::array<bool, 8>) -> uint8_t;
auto unpack_8_booleans(uint8_t) -> std::array<bool, 8>;

// Read from and write to a file
// Note: not using references for `filename` to allow a c-style string literal to be passed in.
auto write_n_bytes(std::string filename, size_t n_bytes, const void* buffer) -> RTNType;
template <typename T>
auto read_whole_file(std::string filename) -> std::vector<T>;

// Upon success, it returns a vector of size `n_bytes`.
// Otherwise, it returns an empty vector.
auto read_n_bytes(std::string filename, size_t n_bytes) -> std::vector<uint8_t>;

// Calculate a suite of statistics.
// Note that arr1 is considered as the ground truth array, so it's the range of
//   arr1 that is used internally for psnr calculations.
// If `omp_nthreads` is zero, then it will use the maximum number of threads.
// The return array contains statistics in the following order:
// ret[0] : RMSE
// ret[1] : L-Infinity
// ret[2] : PSNR
// ret[3] : min of arr1
// ret[4] : max of arr1
template <typename T>
auto calc_stats(const T* arr1, const T* arr2, size_t arr_len, size_t omp_nthreads = 0)
    -> std::array<T, 5>;

template <typename T>
auto kahan_summation(const T*, size_t) -> T;

// Given a whole volume size and a desired chunk size, this helper function
// returns a list of chunks specified by 6 integers:
// chunk[0], [2], [4]: starting index of this chunk in X, Y, and Z;
// chunk[1], [3], [5]: length of this chunk in X, Y, and Z.
// Note 1: the values in `chunk_dim` is only suggestive, meaning that when the volume
//         dimension is not exact multiplies of requested chunk dimension,
//         approximate values are used.
// Note 2: this function works on degraded 2D or 1D volumes too.
auto chunk_volume(const dims_type& vol_dim, const dims_type& chunk_dim)
    -> std::vector<std::array<size_t, 6>>;

// Gather a chunk from a bigger volume
// If the requested chunk lives outside of the volume, whole or part,
// this function returns an empty vector.
template <typename T1, typename T2>
auto gather_chunk(const T1* vol, dims_type vol_dim, const std::array<size_t, 6>& chunk)
    -> std::vector<T2>;

// Put this chunk to a bigger volume
// The `big_vol` should have enough space allocated, and the `small_vol` should contain
// enough elements to scatter. Memory errors will occur if the conditions are not met.
template <typename TBIG, typename TSML>
void scatter_chunk(std::vector<TBIG>& big_vol,
                   dims_type vol_dim,
                   const std::vector<TSML>& small_vol,
                   const std::array<size_t, 6>& chunk);

// Structure that holds information extracted from SPERR headers.
// This structure is returned by helper function `parse_header()`.
struct HeaderInfo {
  uint8_t version_major = 0;
  bool zstd_applied = false;
  bool is_3d = false;
  bool orig_is_float = false;

  // This is the dimension of a 3D volume (NX, NY, NZ) or a 2D slice (NX, NY, 1).
  dims_type vol_dims = {0, 0, 0};

  // If the bitstream represents a 3D volume, this field holds the dimension of chunks.
  // For 2D slices, this field holds undefined values.
  dims_type chunk_dims = {0, 0, 0};
};
auto parse_header(const void*) -> HeaderInfo;

// Calculate the variance of a given array.
// In case of arrays of size zero, it will return std::numeric_limits<T>::infinity().
// In case of `omp_nthreads == 0`, it will use all available OpenMP threads.
template <typename T>
auto calc_variance(const T*, size_t len, size_t omp_nthreads = 0) -> T;

// Decide compression mode based on a collection of parameters.
auto compression_mode(size_t bit_budget, double psnr, double pwe) -> CompMode;

};  // namespace sperr


auto sperr::num_of_xforms(size_t len) -> size_t
{
  assert(len > 0);
  // I decide 8 is the minimal length to do one level of xform.
  const auto f = std::log2(double(len) / 8.0);
  const auto num = f < 0.0 ? size_t{0} : static_cast<size_t>(f) + 1;
  // I also decide that no matter what the input size is,
  // six (6) is the maxinum number of transforms to do.
  return std::min(num, size_t{6});
}

auto sperr::num_of_partitions(size_t len) -> size_t
{
  size_t num_of_parts = 0;  // Num. of partitions we can do
  while (len > 1) {
    num_of_parts++;
    len -= len / 2;
  }

  return num_of_parts;
}

auto sperr::calc_approx_detail_len(size_t orig_len, size_t lev) -> std::array<size_t, 2>
{
  size_t low_len = orig_len;
  size_t high_len = 0;
  for (size_t i = 0; i < lev; i++) {
    high_len = low_len / 2;
    low_len -= high_len;
  }

  return {low_len, high_len};
}

// Good solution to deal with bools and unsigned chars
// https://stackoverflow.com/questions/8461126/how-to-create-a-byte-out-of-8-bool-values-and-vice-versa
auto sperr::pack_booleans(std::vector<uint8_t>& dest, const std::vector<bool>& src, size_t offset)
    -> RTNType
{
  // `src` has to have a size of multiples of 8.
  if (src.size() % 8 != 0)
    return RTNType::BitstreamWrongLen;

  // `dest` should have enough space, as the API specifies.
  assert(dest.size() >= offset + src.size() / 8);

  // How many bits to process at a time.
  constexpr uint64_t bit_stride = 2048;
  constexpr uint64_t byte_stride = bit_stride / 8;
  constexpr uint64_t magic = 0x8040201008040201;

  // It turns out C++ doesn't specify bool to be 1 byte, so we have to use
  // uint8_t here which is definitely 1 byte in size.
  // Also, C++ guarantees conversion between booleans and integers:
  // true <--> 1, false <--> 0.
  auto a = std::array<uint8_t, bit_stride>();
  auto t = std::array<uint64_t, byte_stride>();
  size_t dest_idx = offset;
  auto itr_finish = src.cbegin() + (src.size() / bit_stride) * bit_stride;

  for (auto itr = src.cbegin(); itr != itr_finish; itr += bit_stride) {
    std::copy(itr, itr + bit_stride, a.begin());
    std::memcpy(t.data(), a.data(), a.size());
    std::transform(t.cbegin(), t.cend(), dest.begin() + dest_idx,
                   [](auto e) { return (magic * e) >> 56; });
    dest_idx += byte_stride;
  }

  // For the remaining bits, process 8 bits at a time.
  for (auto itr = itr_finish; itr != src.cend(); itr += 8) {
    std::copy(itr, itr + 8, a.begin());
    std::memcpy(t.data(), a.data(), 8);
    dest[dest_idx++] = (magic * t[0]) >> 56;
  }

  return RTNType::Good;
}

auto sperr::unpack_booleans(std::vector<bool>& dest,
                            const void* src,
                            size_t src_len,
                            size_t src_offset) -> RTNType
{
  // For some reason, a strided unpack implementation similar to the striding
  // strategy used in `pack_booleans()` would result in slower runtime...
  // See Github issue #122.

  if (src == nullptr)
    return RTNType::InvalidParam;

  if (src_len < src_offset)
    return RTNType::BitstreamWrongLen;

  const size_t num_of_bytes = src_len - src_offset;

  // `dest` needs to have enough space allocated, as the API specifies.
  assert(dest.size() >= num_of_bytes * 8);

  const uint8_t* src_ptr = static_cast<const uint8_t*>(src) + src_offset;
  const uint64_t magic = 0x8040201008040201;
  const uint64_t mask = 0x8080808080808080;

#ifndef OMP_UNPACK_BOOLEANS
  // Serial implementation
  //
  auto a = std::array<uint8_t, 8>();
  auto t = uint64_t{0};
  auto b8 = uint64_t{0};
  auto dest_itr = dest.begin();
  for (size_t byte_idx = 0; byte_idx < num_of_bytes; byte_idx++) {
    b8 = *(src_ptr + byte_idx);
    t = ((magic * b8) & mask) >> 7;
    std::memcpy(a.data(), &t, 8);
    std::copy(a.cbegin(), a.cend(), dest_itr);
    dest_itr += 8;
  }
#else
  // Parallel implementation
  //
  // Because in most implementations std::vector<bool> is stored as uint64_t
  // values, we parallel in strides of 64 bits, or 8 bytes.
  const size_t stride_size = 8;
  const size_t num_of_strides = num_of_bytes / stride_size;

#pragma omp parallel for
  for (size_t stride = 0; stride < num_of_strides; stride++) {
    uint8_t a[64]{0};
    for (size_t byte = 0; byte < stride_size; byte++) {
      const uint8_t* ptr = src_ptr + stride * stride_size + byte;
      const uint64_t t = ((magic * (*ptr)) & mask) >> 7;
      std::memcpy(a + byte * 8, &t, 8);
    }
    for (size_t i = 0; i < 64; i++)
      dest[stride * 64 + i] = a[i];
  }
  // This loop is at most 7 iterations, so not to worry about parallel anymore.
  for (size_t byte_idx = stride_size * num_of_strides; byte_idx < num_of_bytes; byte_idx++) {
    const uint8_t* ptr = src_ptr + byte_idx;
    const uint64_t t = ((magic * (*ptr)) & mask) >> 7;
    uint8_t a[8]{0};
    std::memcpy(a, &t, 8);
    for (size_t i = 0; i < 8; i++)
      dest[byte_idx * 8 + i] = a[i];
  }
#endif

  return RTNType::Good;
}

auto sperr::pack_8_booleans(std::array<bool, 8> src) -> uint8_t
{
  // It turns out that C++ doesn't specify bool to be one byte,
  // so to be safe we copy the content of src to array of uint8_t.
  auto bytes = std::array<uint8_t, 8>();
  std::copy(src.begin(), src.end(), bytes.begin());
  const uint64_t magic = 0x8040201008040201;
  uint64_t t = 0;
  std::memcpy(&t, bytes.data(), 8);
  uint8_t dest = (magic * t) >> 56;
  return dest;
}

auto sperr::unpack_8_booleans(uint8_t src) -> std::array<bool, 8>
{
  const uint64_t magic = 0x8040201008040201;
  const uint64_t mask = 0x8080808080808080;
  uint64_t t = ((magic * src) & mask) >> 7;
  // It turns out that C++ doesn't specify bool to be one byte,
  // so to be safe we use an array of uint8_t.
  auto bytes = std::array<uint8_t, 8>();
  std::memcpy(bytes.data(), &t, 8);
  auto b8 = std::array<bool, 8>();
  std::copy(bytes.begin(), bytes.end(), b8.begin());
  return b8;
}

auto sperr::read_n_bytes(std::string filename, size_t n_bytes) -> std::vector<uint8_t>
{
  auto buf = std::vector<uint8_t>();

  std::unique_ptr<std::FILE, decltype(&std::fclose)> fp(std::fopen(filename.data(), "rb"),
                                                        &std::fclose);

  if (!fp)
    return buf;

  std::fseek(fp.get(), 0, SEEK_END);
  if (std::ftell(fp.get()) < n_bytes)
    return buf;

  std::fseek(fp.get(), 0, SEEK_SET);
  buf.resize(n_bytes);
  if (std::fread(buf.data(), 1, n_bytes, fp.get()) != n_bytes)
    buf.clear();

  return buf;
}

template <typename T>
auto sperr::read_whole_file(std::string filename) -> std::vector<T>
{
  std::vector<T> buf;

  std::unique_ptr<std::FILE, decltype(&std::fclose)> fp(std::fopen(filename.data(), "rb"),
                                                        &std::fclose);
  if (!fp)
    return buf;

  std::fseek(fp.get(), 0, SEEK_END);
  const size_t file_size = std::ftell(fp.get());
  const size_t num_vals = file_size / sizeof(T);
  std::fseek(fp.get(), 0, SEEK_SET);

  buf.resize(num_vals);
  size_t nread = std::fread(buf.data(), sizeof(T), num_vals, fp.get());
  if (nread != num_vals)
    buf.clear();

  return buf;
}
template auto sperr::read_whole_file(std::string) -> std::vector<float>;
template auto sperr::read_whole_file(std::string) -> std::vector<double>;
template auto sperr::read_whole_file(std::string) -> std::vector<uint8_t>;

auto sperr::write_n_bytes(std::string filename, size_t n_bytes, const void* buffer) -> RTNType
{
  std::unique_ptr<std::FILE, decltype(&std::fclose)> fp(std::fopen(filename.data(), "wb"),
                                                        &std::fclose);
  if (!fp)
    return RTNType::IOError;

  if (std::fwrite(buffer, 1, n_bytes, fp.get()) != n_bytes)
    return RTNType::IOError;
  else
    return RTNType::Good;
}

template <typename T>
auto sperr::calc_stats(const T* arr1, const T* arr2, size_t arr_len, size_t omp_nthreads)
    -> std::array<T, 5>
{
  const size_t stride_size = 4096;
  const size_t num_of_strides = arr_len / stride_size;
  const size_t remainder_size = arr_len - stride_size * num_of_strides;

  // Use the maximum possible threads if 0 is passed in.
#ifdef USE_OMP
  if (omp_nthreads == 0)
    omp_nthreads = omp_get_max_threads();
#endif

  auto rmse = T{0.0};
  auto linfty = T{0.0};
  auto psnr = T{0.0};
  auto arr1min = T{0.0};
  auto arr1max = T{0.0};

  //
  // Calculate min and max of arr1
  //
  const auto minmax = std::minmax_element(arr1, arr1 + arr_len);
  arr1min = *minmax.first;
  arr1max = *minmax.second;

  //
  // In rare cases, the two input arrays are identical.
  //
  auto is_equal = std::equal(arr1, arr1 + arr_len, arr2);
  if (is_equal) {
    psnr = std::numeric_limits<T>::infinity();
    return {rmse, linfty, psnr, arr1min, arr1max};
  }

  auto sum_vec = std::vector<T>(num_of_strides + 1);
  auto linfty_vec = std::vector<T>(num_of_strides + 1);

//
// Calculate diff summation and l-infty of each stride
//
#pragma omp parallel for num_threads(omp_nthreads)
  for (size_t stride_i = 0; stride_i < num_of_strides; stride_i++) {
    T maxerr = 0.0;
    auto buf = std::array<T, stride_size>();
    for (size_t i = 0; i < stride_size; i++) {
      const size_t idx = stride_i * stride_size + i;
      auto diff = std::abs(arr1[idx] - arr2[idx]);
      maxerr = std::max(maxerr, diff);
      buf[i] = diff * diff;
    }
    sum_vec[stride_i] = std::accumulate(buf.cbegin(), buf.cend(), T{0.0});
    linfty_vec[stride_i] = maxerr;
  }

  //
  // Calculate diff summation and l-infty of the remaining elements
  //
  T last_linfty = 0.0;
  auto last_buf = std::array<T, stride_size>{};  // must be enough for `remainder_size` elements.
  for (size_t i = 0; i < remainder_size; i++) {
    const size_t idx = stride_size * num_of_strides + i;
    auto diff = std::abs(arr1[idx] - arr2[idx]);
    last_linfty = std::max(last_linfty, diff);
    last_buf[i] = diff * diff;
  }
  sum_vec[num_of_strides] = T{0.0};
  sum_vec[num_of_strides] =
      std::accumulate(last_buf.cbegin(), last_buf.cbegin() + remainder_size, T{0.0});
  linfty_vec[num_of_strides] = last_linfty;

  //
  // Now calculate linfty
  //
  linfty = *(std::max_element(linfty_vec.begin(), linfty_vec.end()));

  //
  // Now calculate rmse and psnr
  // Note: psnr is calculated in dB, and follows the equation described in:
  // http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/VELDHUIZEN/node18.html
  // Also refer to https://www.mathworks.com/help/vision/ref/psnr.html
  //
  const auto mse = std::accumulate(sum_vec.cbegin(), sum_vec.cend(), T{0.0}) / T(arr_len);
  rmse = std::sqrt(mse);
  const auto range_sq = (arr1max - arr1min) * (arr1max - arr1min);
  psnr = std::log10(range_sq / mse) * T{10.0};

  return {rmse, linfty, psnr, arr1min, arr1max};
}
template auto sperr::calc_stats(const float*, const float*, size_t, size_t) -> std::array<float, 5>;
template auto sperr::calc_stats(const double*, const double*, size_t, size_t)
    -> std::array<double, 5>;

template <typename T>
auto sperr::kahan_summation(const T* arr, size_t len) -> T
{
  T sum = 0.0, c = 0.0;
  T t, y;
  for (size_t i = 0; i < len; i++) {
    y = arr[i] - c;
    t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }

  return sum;
}
template auto sperr::kahan_summation(const float*, size_t) -> float;
template auto sperr::kahan_summation(const double*, size_t) -> double;

auto sperr::chunk_volume(const std::array<size_t, 3>& vol_dim,
                         const std::array<size_t, 3>& chunk_dim)
    -> std::vector<std::array<size_t, 6>>
{
  // Step 1: figure out how many segments are there along each axis.
  auto n_segs = std::array<size_t, 3>();
  for (size_t i = 0; i < 3; i++) {
    n_segs[i] = vol_dim[i] / chunk_dim[i];
    // In case the last segment is shorter than `chunk_dim[i]`, test if it's
    // longer than half of `chunk_dim[i]`. If it is, it can qualify another segment.
    if ((vol_dim[i] % chunk_dim[i]) > (chunk_dim[i] / 2))
      n_segs[i]++;
    // In case the volume has one dimension that's too small, let's make it
    // at least one segment anyway.
    if (n_segs[i] == 0)
      n_segs[i] = 1;
  }

  // Step 2: calculate the starting indices of each segment along each axis.
  auto x_tics = std::vector<size_t>(n_segs[0] + 1);
  for (size_t i = 0; i < n_segs[0]; i++)
    x_tics[i] = i * chunk_dim[0];
  x_tics[n_segs[0]] = vol_dim[0];  // the last tic is the length in X

  auto y_tics = std::vector<size_t>(n_segs[1] + 1);
  for (size_t i = 0; i < n_segs[1]; i++)
    y_tics[i] = i * chunk_dim[1];
  y_tics[n_segs[1]] = vol_dim[1];  // the last tic is the length in Y

  auto z_tics = std::vector<size_t>(n_segs[2] + 1);
  for (size_t i = 0; i < n_segs[2]; i++)
    z_tics[i] = i * chunk_dim[2];
  z_tics[n_segs[2]] = vol_dim[2];  // the last tic is the length in Z

  // Step 3: fill in details of each chunk
  auto n_chunks = n_segs[0] * n_segs[1] * n_segs[2];
  auto chunks = std::vector<std::array<size_t, 6>>(n_chunks);
  size_t chunk_idx = 0;
  for (size_t z = 0; z < n_segs[2]; z++)
    for (size_t y = 0; y < n_segs[1]; y++)
      for (size_t x = 0; x < n_segs[0]; x++) {
        chunks[chunk_idx][0] = x_tics[x];                  // X start
        chunks[chunk_idx][1] = x_tics[x + 1] - x_tics[x];  // X length
        chunks[chunk_idx][2] = y_tics[y];                  // Y start
        chunks[chunk_idx][3] = y_tics[y + 1] - y_tics[y];  // Y length
        chunks[chunk_idx][4] = z_tics[z];                  // Z start
        chunks[chunk_idx][5] = z_tics[z + 1] - z_tics[z];  // Z length
        chunk_idx++;
      }

  return chunks;
}

template <typename T1, typename T2>
auto sperr::gather_chunk(const T1* vol, dims_type vol_dim, const std::array<size_t, 6>& chunk)
    -> std::vector<T2>
{
  auto chunk_buf = std::vector<T2>();
  if (chunk[0] + chunk[1] > vol_dim[0] || chunk[2] + chunk[3] > vol_dim[1] ||
      chunk[4] + chunk[5] > vol_dim[2])
    return chunk_buf;

  chunk_buf.resize(chunk[1] * chunk[3] * chunk[5]);

  size_t idx = 0;
  for (size_t z = chunk[4]; z < chunk[4] + chunk[5]; z++) {
    const size_t plane_offset = z * vol_dim[0] * vol_dim[1];
    for (size_t y = chunk[2]; y < chunk[2] + chunk[3]; y++) {
      const size_t col_offset = plane_offset + y * vol_dim[0];
      for (size_t x = chunk[0]; x < chunk[0] + chunk[1]; x++)
        chunk_buf[idx++] = vol[col_offset + x];
    }
  }

  // Will be subject to Named Return Value Optimization.
  return chunk_buf;
}
template auto sperr::gather_chunk(const float*, dims_type, const std::array<size_t, 6>&)
    -> std::vector<float>;
template auto sperr::gather_chunk(const float*, dims_type, const std::array<size_t, 6>&)
    -> std::vector<double>;
template auto sperr::gather_chunk(const double*, dims_type, const std::array<size_t, 6>&)
    -> std::vector<float>;
template auto sperr::gather_chunk(const double*, dims_type, const std::array<size_t, 6>&)
    -> std::vector<double>;

template <typename TBIG, typename TSML>
void sperr::scatter_chunk(std::vector<TBIG>& big_vol,
                          dims_type vol_dim,
                          const std::vector<TSML>& small_vol,
                          const std::array<size_t, 6>& chunk)
{
  size_t idx = 0;
  for (size_t z = chunk[4]; z < chunk[4] + chunk[5]; z++) {
    const size_t plane_offset = z * vol_dim[0] * vol_dim[1];
    for (size_t y = chunk[2]; y < chunk[2] + chunk[3]; y++) {
      const size_t col_offset = plane_offset + y * vol_dim[0];
      for (size_t x = chunk[0]; x < chunk[0] + chunk[1]; x++)
        big_vol[col_offset + x] = small_vol[idx++];
    }
  }
}
template void sperr::scatter_chunk(std::vector<float>&,
                                   dims_type,
                                   const std::vector<float>&,
                                   const std::array<size_t, 6>&);
template void sperr::scatter_chunk(std::vector<float>&,
                                   dims_type,
                                   const std::vector<double>&,
                                   const std::array<size_t, 6>&);
template void sperr::scatter_chunk(std::vector<double>&,
                                   dims_type,
                                   const std::vector<float>&,
                                   const std::array<size_t, 6>&);
template void sperr::scatter_chunk(std::vector<double>&,
                                   dims_type,
                                   const std::vector<double>&,
                                   const std::array<size_t, 6>&);

auto sperr::parse_header(const void* ptr) -> HeaderInfo
{
  const uint8_t* u8p = static_cast<const uint8_t*>(ptr);
  size_t loc = 0;
  auto header = HeaderInfo();

  // Parse version numbers
  header.version_major = *u8p;
  loc++;

  // Parse 8 booleans
  const auto b8 = sperr::unpack_8_booleans(u8p[loc]);
  loc++;

  header.zstd_applied = b8[0];
  header.is_3d = b8[1];
  header.orig_is_float = b8[2];

  // Parse the dimension info.
  if (header.is_3d) {
    if (b8[3]) {  // there are multiple chunks!
      uint32_t vcdim[6];
      std::memcpy(vcdim, u8p + loc, sizeof(vcdim));
      loc += sizeof(vcdim);

      header.vol_dims[0] = vcdim[0];
      header.vol_dims[1] = vcdim[1];
      header.vol_dims[2] = vcdim[2];
      header.chunk_dims[0] = vcdim[3];
      header.chunk_dims[1] = vcdim[4];
      header.chunk_dims[2] = vcdim[5];
    }
    else {
      uint32_t vdim[3];
      std::memcpy(vdim, u8p + loc, sizeof(vdim));
      loc += sizeof(vdim);

      header.vol_dims[0] = vdim[0];
      header.vol_dims[1] = vdim[1];
      header.vol_dims[2] = vdim[2];
      header.chunk_dims[0] = vdim[0];
      header.chunk_dims[1] = vdim[1];
      header.chunk_dims[2] = vdim[2];
    }
  }
  else {
    uint32_t dims[2];
    std::memcpy(dims, u8p + loc, sizeof(dims));
    loc += sizeof(dims);

    header.vol_dims[0] = dims[0];
    header.vol_dims[1] = dims[1];
    header.vol_dims[2] = 1;
  }

  return header;
}

template <typename T>
auto sperr::calc_variance(const T* arr, size_t len, size_t omp_nthreads) -> T
{
  if (len == 0)
    return std::numeric_limits<T>::infinity();

#ifdef USE_OMP
  if (omp_nthreads == 0)
    omp_nthreads = omp_get_max_threads();
#endif

  const size_t stride_size = 4096;
  const size_t num_strides = len / stride_size;
  auto tmp_buf = std::vector<T>(num_strides + 1);

  // First, calculate the mean of this array.
#pragma omp parallel for num_threads(omp_nthreads)
  for (size_t i = 0; i < num_strides; i++) {
    const T* beg = arr + i * stride_size;
    tmp_buf[i] = std::accumulate(beg, beg + stride_size, T{0.0});
  }
  tmp_buf[num_strides] = 0.0;
  tmp_buf[num_strides] = std::accumulate(arr + num_strides * stride_size, arr + len, T{0.0});
  const auto elem_sum = std::accumulate(tmp_buf.cbegin(), tmp_buf.cend(), T{0.0});
  const auto mean = elem_sum / static_cast<T>(len);

  // Second, calculate the variance.
#pragma omp parallel for num_threads(omp_nthreads)
  for (size_t i = 0; i < num_strides; i++) {
    const T* beg = arr + i * stride_size;
    tmp_buf[i] = std::accumulate(beg, beg + stride_size, T{0.0}, [mean](auto init, auto v) {
      return init + (v - mean) * (v - mean);
    });
  }
  tmp_buf[num_strides] = 0.0;
  tmp_buf[num_strides] =
      std::accumulate(arr + num_strides * stride_size, arr + len, T{0.0},
                      [mean](auto init, auto v) { return init + (v - mean) * (v - mean); });
  const auto diff_sum = std::accumulate(tmp_buf.cbegin(), tmp_buf.cend(), T{0.0});
  const auto var = diff_sum / static_cast<T>(len);

  return var;
}
template auto sperr::calc_variance(const float*, size_t, size_t) -> float;
template auto sperr::calc_variance(const double*, size_t, size_t) -> double;

auto sperr::compression_mode(size_t bit_budget, double psnr, double pwe) -> CompMode
{
  if (bit_budget < sperr::max_size && psnr == sperr::max_d && pwe == 0.0) {
    return CompMode::FixedSize;
  }
  else if (bit_budget == sperr::max_size && psnr < sperr::max_d && pwe == 0.0) {
    return CompMode::FixedPSNR;
  }
  else if (bit_budget == sperr::max_size && psnr == sperr::max_d && pwe > 0.0) {
    return CompMode::FixedPWE;
  }
  else {
    return CompMode::Unknown;
  }
}



#endif
