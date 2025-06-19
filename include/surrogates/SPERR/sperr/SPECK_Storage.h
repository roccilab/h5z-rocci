#ifndef SPECK_STORAGE_H
#define SPECK_STORAGE_H

//
// This class implements mechanisms and interfaces to help 2D and 3D SPECK
// classes manage storage of coefficients and bitstreams.
//


#include "sperr_helper.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <numeric>

using d2_type = std::array<double, 2>;
namespace sperr {

class SPECK_Storage {
 public:
  //
  // Memory management: input and output
  //
  template <typename T>
  auto copy_data(const T*, size_t len, dims_type dims) -> RTNType;
  auto take_data(vecd_type&&, dims_type dims) -> RTNType;
  auto release_data() -> vecd_type&&;
  auto view_data() const -> const vecd_type&;
  auto view_encoded_bitstream() const -> const vec8_type&;
  auto release_encoded_bitstream() -> vec8_type&&;
  auto release_quantized_coeff() -> vecd_type&&;
  size_t count_LSP();

  // Prepare internal states for a decompression operation from an encoded bitstream
  // Note: it takes a raw pointer because it accesses memory provided by others,
  //       and others most likely provide a raw pointer.
  auto parse_encoded_bitstream(const void*, size_t) -> RTNType;

  // Given a SPECK stream, tell how long the speck stream (including header) is.
  // Note: don't need to provide the buffer size because this function
  //       goes to a fixed location to retrieve the stream size.
  auto get_speck_stream_size(const void*) const -> uint64_t;

  void set_dimensions(dims_type);
  void set_data_range(double);

  // Set all compression parameters together. Notice their order!
  auto set_comp_params(size_t bit_budget, double psnr, double pwe) -> RTNType;

 protected:
  //
  // Member variables
  //
  const size_t m_header_size = 16;  // See header definition in SPECK_Storage.cpp.
  const size_t m_u64_garbage_val = std::numeric_limits<size_t>::max();
  size_t m_encode_budget = 0;
  size_t m_LSP_mask_cnt = 0;           // Number of TRUE values in `m_LSP_mask`. Decoding only
  size_t m_bit_idx = 0;                // Which bit we're at? Decoding only
  double m_max_threshold = 0.0;       // float representation of max threshold
  double m_data_range = sperr::max_d;  // range of data before DWT
  double m_target_pwe = 0.0;           // used in fixed-PWE mode
  double m_threshold = 0.0;            // Threshold that's used for an iteration
  double m_target_psnr = sperr::max_d;

  // Use a cache variable to store the current compression mode.
  // The compression mode, however, is NOT set or determined by `m_mode_cache`. It is determined
  // by the combination of `m_encode_budget`, `m_qz_lev`, `m_target_psnr`, and `m_target_pwe`,
  // and can be retrieved by calling function sperr::compression_mode().
  CompMode m_mode_cache = CompMode::Unknown;

  // Data storage shared by both 2D and 3D SPECK.
  //
  std::vector<double> m_coeff_buf;  // Wavelet coefficients
  std::vector<bool> m_bit_buffer;   // Bitstream produced by the algorithm
  std::vector<size_t> m_LIP;        // List of insignificant pixels
  std::vector<size_t> m_LSP_new;    // List of newly found significant pixels
  std::vector<bool> m_LSP_mask;     // Significant pixels previously found
  std::vector<bool> m_sign_array;   // Keep the signs of every coefficient
  vec8_type m_encoded_stream;       // Stores the SPECK bitstream
  dims_type m_dims = {0, 0, 0};     // Dimension of the 2D/3D volume

  // In fixed-pwe mode, keep track of would-be quantized coefficient values.
  // These values are reconstructed exactly the same way as if during decompression.
  // Other mathematically equivalent methods would however not give the same values.
  std::vector<double> m_qz_coeff;

  //
  // Member methods
  //
  auto m_prepare_encoded_bitstream() -> RTNType;
  auto m_refinement_pass_encode() -> RTNType;
  auto m_refinement_pass_decode() -> RTNType;

  auto m_estimate_rmse(double q) const -> double;
  auto m_estimate_finest_q(const double &eb_coeff) const -> double;
};

};  // namespace sperr


template <typename T>
auto sperr::SPECK_Storage::copy_data(const T* p, size_t len, dims_type dims) -> RTNType
{
  static_assert(std::is_floating_point<T>::value, "!! Only floating point values are supported !!");
  if (len != dims[0] * dims[1] * dims[2])
    return RTNType::WrongDims;

  m_coeff_buf.resize(len);
  std::copy(p, p + len, m_coeff_buf.begin());

  m_dims = dims;

  return RTNType::Good;
}
template auto sperr::SPECK_Storage::copy_data(const double*, size_t, dims_type) -> RTNType;
template auto sperr::SPECK_Storage::copy_data(const float*, size_t, dims_type) -> RTNType;

auto sperr::SPECK_Storage::take_data(vecd_type&& coeffs, dims_type dims) -> RTNType
{
  if (coeffs.size() != dims[0] * dims[1] * dims[2])
    return RTNType::WrongDims;

  m_coeff_buf = std::move(coeffs);
  m_dims = dims;

  return RTNType::Good;
}

auto sperr::SPECK_Storage::release_data() -> vecd_type&&
{
  m_dims = {0, 0, 0};
  return std::move(m_coeff_buf);
}

auto sperr::SPECK_Storage::view_data() const -> const vecd_type&
{
  return m_coeff_buf;
}

auto sperr::SPECK_Storage::release_quantized_coeff() -> vecd_type&&
{
  if (!m_qz_coeff.empty()) {
    assert(m_qz_coeff.size() == m_sign_array.size());
    const auto tmp = d2_type{-1.0, 1.0};
    std::transform(m_qz_coeff.begin(), m_qz_coeff.end(), m_sign_array.cbegin(), m_qz_coeff.begin(),
                   [tmp](auto v, auto b) { return v * tmp[b]; });
  }

  return std::move(m_qz_coeff);
}

void sperr::SPECK_Storage::set_data_range(double range)
{
  m_data_range = range;
}

auto sperr::SPECK_Storage::m_prepare_encoded_bitstream() -> RTNType
{
  // Header definition: 16 bytes in total:
  // m_max_threshold,   num_useful_bits
  // double,             uint64_t

  const uint64_t useful_bits = m_bit_buffer.size();
  while (m_bit_buffer.size() % 8 != 0)
    m_bit_buffer.push_back(false);
  const size_t total_size = m_header_size + m_bit_buffer.size() / 8;
  m_encoded_stream.resize(total_size);
  

  // Fill header
  size_t pos = 0;
  auto* const ptr = m_encoded_stream.data();
  std::memcpy(ptr + pos, &m_max_threshold, sizeof(m_max_threshold));
  pos += sizeof(m_max_threshold);

  std::memcpy(ptr + pos, &useful_bits, sizeof(useful_bits));
  pos += sizeof(useful_bits);
  assert(pos == m_header_size);

  // Assemble the bitstream into bytes
  auto rtn = sperr::pack_booleans(m_encoded_stream, m_bit_buffer, pos);
  // Restore the number of useful bits
  m_bit_buffer.resize(useful_bits);
  return rtn;
}

auto sperr::SPECK_Storage::parse_encoded_bitstream(const void* comp_buf, size_t comp_size)
    -> RTNType
{
  // The buffer passed in is supposed to consist a header and then a compacted
  // bitstream, just like what was returned by `prepare_encoded_bitstream()`. Note:
  // header definition is documented in prepare_encoded_bitstream().

  const uint8_t* const ptr = static_cast<const uint8_t*>(comp_buf);

  // Parse the header
  size_t pos = 0;
  std::memcpy(&m_max_threshold, ptr + pos, sizeof(m_max_threshold));
  pos += sizeof(m_max_threshold);

  uint64_t useful_bits = 0;
  std::memcpy(&useful_bits, ptr + pos, sizeof(useful_bits));
  pos += sizeof(useful_bits);
  auto num_bytes = useful_bits / 8;
  if (useful_bits % 8 != 0)
    ++num_bytes;
  if (num_bytes != comp_size - pos)
    return RTNType::BitstreamWrongLen;

  // Unpack bits
  m_bit_buffer.resize(num_bytes * 8);  // Allocate enough space before passing it around
  auto rtn = sperr::unpack_booleans(m_bit_buffer, comp_buf, comp_size, pos);
  m_bit_buffer.resize(useful_bits);  // Restore the actual useful number of bits

  return rtn;
}

auto sperr::SPECK_Storage::view_encoded_bitstream() const -> const vec8_type&
{
  return m_encoded_stream;
}

auto sperr::SPECK_Storage::release_encoded_bitstream() -> vec8_type&&
{
  return std::move(m_encoded_stream);
}

auto sperr::SPECK_Storage::get_speck_stream_size(const void* buf) const -> uint64_t
{
  // Given the header definition in `prepare_encoded_bitstream()`, directly
  // go retrieve the value stored in the last 8 bytes of the header
  //
  const uint8_t* const ptr = static_cast<const uint8_t*>(buf);
  uint64_t useful_bits;
  std::memcpy(&useful_bits, ptr + m_header_size - 8, sizeof(useful_bits));

  auto num_bytes = useful_bits / 8;
  if (useful_bits % 8 != 0)
    ++num_bytes;

  return (m_header_size + num_bytes);
}

void sperr::SPECK_Storage::set_dimensions(dims_type dims)
{
  m_dims = dims;
}

auto sperr::SPECK_Storage::set_comp_params(size_t budget, double psnr, double pwe) -> RTNType
{
  // First set those ones that only need a plain copy
  m_target_psnr = psnr;
  m_target_pwe = pwe;

  // Second set bit budget, which would require a little manipulation.
  if (budget == sperr::max_size) {
    m_encode_budget = sperr::max_size;
    return RTNType::Good;
  }

  if (budget <= m_header_size * 8) {
    m_encode_budget = 0;
    return RTNType::Error;
  }
  budget -= m_header_size * 8;

  assert(budget % 8 == 0);
  m_encode_budget = budget;

  return RTNType::Good;
}

auto sperr::SPECK_Storage::m_refinement_pass_encode() -> RTNType
{
  // First, process significant pixels previously found.
  //
  const auto tmpd = d2_type{m_threshold * -0.5, m_threshold * 0.5};
  const auto tmptr = d2_type{0.0, -m_threshold};

  assert(m_encode_budget >= m_bit_buffer.size());

  if (m_mode_cache == CompMode::FixedPWE) {
    for (size_t i = 0; i < m_LSP_mask.size(); i++) {
      if (m_LSP_mask[i]) {
        const bool o1 = m_coeff_buf[i] >= m_threshold;
        m_bit_buffer.push_back(o1);
        m_coeff_buf[i] += tmptr[o1];
        m_qz_coeff[i] += tmpd[o1];
      }
    }
  }
  else {
    for (size_t i = 0; i < m_LSP_mask.size(); i++) {
      if (m_LSP_mask[i]) {
        const bool o1 = m_coeff_buf[i] >= m_threshold;
        m_bit_buffer.push_back(o1);
        m_coeff_buf[i] += tmptr[o1];
      }
    }
  }

  // If we generated more than necessary bits (in fixed-rate mode), then we resize!
  if (m_bit_buffer.size() >= m_encode_budget) {
    m_bit_buffer.resize(m_encode_budget);
    return RTNType::BitBudgetMet;
  }

  // Second, mark newly found significant pixels in `m_LSP_mask`.
  //
  for (auto idx : m_LSP_new)
    m_LSP_mask[idx] = true;
  m_LSP_new.clear();

  return RTNType::Good;
}

auto sperr::SPECK_Storage::m_refinement_pass_decode() -> RTNType
{
  // First, process significant pixels previously found.
  //
  const auto tmpd = d2_type{m_threshold * -0.5, m_threshold * 0.5};

  assert(m_bit_buffer.size() >= m_bit_idx);
  if (m_bit_buffer.size() - m_bit_idx > m_LSP_mask_cnt) {  // No need to check BitBudgetMet
    for (size_t i = 0; i < m_LSP_mask.size(); i++) {
      if (m_LSP_mask[i])
        m_coeff_buf[i] += tmpd[m_bit_buffer[m_bit_idx++]];
    }
  }
  else {  // Need to check BitBudgetMet
    for (size_t i = 0; i < m_LSP_mask.size(); i++) {
      if (m_LSP_mask[i]) {
        m_coeff_buf[i] += tmpd[m_bit_buffer[m_bit_idx++]];
        if (m_bit_idx >= m_bit_buffer.size())
          return RTNType::BitBudgetMet;
      }
    }
  }

  // Second, mark newly found significant pixels in `m_LSP_mark`
  //
  for (auto idx : m_LSP_new)
    m_LSP_mask[idx] = true;
  m_LSP_mask_cnt += m_LSP_new.size();
  m_LSP_new.clear();

  return RTNType::Good;
}

auto sperr::SPECK_Storage::m_estimate_rmse(double q) const -> double
{
  const auto half_q = q * 0.5;
  const auto len = m_coeff_buf.size();
  const size_t stride_size = 4096;
  const size_t num_strides = len / stride_size;
  //const size_t remainder_size = len - stride_size * num_strides;
  auto tmp_buf = vecd_type(num_strides + 1);

  for (size_t i = 0; i < num_strides; i++) {
    const auto beg = m_coeff_buf.cbegin() + i * stride_size;
    tmp_buf[i] = std::accumulate(beg, beg + stride_size, 0.0, [half_q, q](auto init, auto v) {
      auto diff = v >= q ? std::remainder(v + half_q, q) : v;
      return init + diff * diff;
    });
  }

  // Let's also process the last stride.
  tmp_buf[num_strides] = 0.0;
  tmp_buf[num_strides] = std::accumulate(m_coeff_buf.cbegin() + num_strides * stride_size,
                                         m_coeff_buf.cend(), 0.0, [half_q, q](auto init, auto v) {
                                           auto diff = v >= q ? std::remainder(v + half_q, q) : v;
                                           return init + diff * diff;
                                         });
  const auto total_sum = std::accumulate(tmp_buf.cbegin(), tmp_buf.cend(), 0.0);
  const auto mse = total_sum / static_cast<double>(len);

  return std::sqrt(mse);
}



auto sperr::SPECK_Storage::m_estimate_finest_q(const double &coeff=1.5) const -> double
{
  if (m_mode_cache == CompMode::FixedPWE) {
    // Most recent tests show that q = 1.5 t is a pretty good estimation.
    return m_target_pwe * coeff;
  }
  else if (m_mode_cache == CompMode::FixedPSNR) {
    // Note: based on Peter's estimation method, to achieved the target PSNR, the terminal
    //       quantization threshold should be (2.0 * sqrt(3.0) * rmse).
    const auto t_mse = (m_data_range * m_data_range) * std::pow(10.0, -m_target_psnr / 10.0);
    const auto t_rmse = std::sqrt(t_mse);
    auto q = 2.0 * std::sqrt(t_mse * 3.0);
    while (m_estimate_rmse(q) > t_rmse)
      q /= std::pow(2.0, 0.25);  // Four adjustments would effectively halve q.
    return q;
  }
  else
    return 0.0;
}

size_t sperr::SPECK_Storage::count_LSP(){
  size_t count=0;
  for(size_t i=0;i<m_LSP_mask.size();i++){
    if(m_LSP_mask[i])
      count++;
  }
  return count;
}

#endif
