#ifndef CONDITIONER_H
#define CONDITIONER_H

//
// Applies conditioning operations to an array of data.
//

#include "sperr_helper.h"

#include <tuple>
#include <algorithm>  // std::all_of()
#include <cassert>
#include <cmath>    // std::sqrt()
#include <cstring>  // std::memcpy()
#include <numeric>  // std::accumulate()
#include <vector>

namespace sperr {

class Conditioner {
 public:
  // Define a few data types
  using meta_type = std::array<uint8_t, 17>;
  using settings_type = std::array<bool, 4>;

  // Booleans recording what operations are applied:
  // bool[0] : subtract mean
  // bool[1] : divide by rms
  // bool[2-3] : unused
  void toggle_all_settings(settings_type);

  // The 17 bytes returned by `condition()`: 1 byte (8 booleans) followed by two
  // doubles. The byte of booleans records what operations are applied:
  // Accordingly, `inverse_condition()` takes a buffer of 17 bytes.
  auto condition(vecd_type&) -> std::pair<RTNType, meta_type>;
  auto inverse_condition(vecd_type&, const meta_type&) -> RTNType;

  // Given an array, test if it's a constant field. The test result (true/false)
  // is returned as boolean and a properly packed meta data if true.
  // Similarly, `parse_constant` takes in a meta data block and returns if it
  // represents a constant field. If true, it also returns the constant value.
  auto test_constant(const vecd_type&) const -> std::pair<bool, meta_type>;
  auto parse_constant(const meta_type&) const -> std::tuple<bool, double, uint64_t>;

 private:
  //
  // what settings are applied?
  //
  std::array<bool, 8> m_settings = {true,    // subtract mean
                                    false,   // divide by rms
                                    false,   // added: skip wavelet?
                                    false,   // unused
                                    false,   // [4]: is this a constant field?
                                    false,   // unused
                                    false,   // unused
                                    false};  // unused

  // Calculations are carried out by strides, which
  // should be a divisor of the input data size.
  size_t m_num_strides = 2048;  // should be good enough for most applications.
  std::vector<double> m_stride_buf;

  // Action items
  // Buffers passed in here are guaranteed to have correct lengths and conditions.
  auto m_calc_mean(const vecd_type& buf) -> double;
  auto m_calc_rms(const vecd_type& buf) -> double;
  // Adjust the value of `m_num_strides` so it'll be a divisor of `len`.
  void m_adjust_strides(size_t len);
};

};  // namespace sperr


void sperr::Conditioner::toggle_all_settings(std::array<bool, 4> b4)
{
  m_settings[0] = b4[0];
  m_settings[1] = b4[1];
}

auto sperr::Conditioner::m_calc_mean(const vecd_type& buf) -> double
{
  assert(buf.size() % m_num_strides == 0);

  m_stride_buf.resize(m_num_strides);
  const size_t stride_size = buf.size() / m_num_strides;

  for (size_t s = 0; s < m_num_strides; s++) {
    auto begin = buf.begin() + stride_size * s;
    auto end = begin + stride_size;
    m_stride_buf[s] = std::accumulate(begin, end, double{0.0}) / static_cast<double>(stride_size);
  }

  double sum = std::accumulate(m_stride_buf.begin(), m_stride_buf.end(), double{0.0});

  return (sum / static_cast<double>(m_stride_buf.size()));
}

auto sperr::Conditioner::m_calc_rms(const vecd_type& buf) -> double
{
  assert(buf.size() % m_num_strides == 0);

  m_stride_buf.resize(m_num_strides);
  const size_t stride_size = buf.size() / m_num_strides;

  for (size_t s = 0; s < m_num_strides; s++) {
    auto begin = buf.begin() + stride_size * s;
    auto end = begin + stride_size;
    m_stride_buf[s] =
        std::accumulate(begin, end, double{0.0}, [](auto a, auto b) { return a + b * b; });

    m_stride_buf[s] /= static_cast<double>(stride_size);
  }

  double rms = std::accumulate(m_stride_buf.begin(), m_stride_buf.end(), double{0.0});
  rms /= static_cast<double>(m_stride_buf.size());
  rms = std::sqrt(rms);

  return rms;
}

void sperr::Conditioner::m_adjust_strides(size_t len)
{
  // Let's say 2048 is a starting point
  m_num_strides = 2048;
  if (len % m_num_strides == 0)
    return;

  size_t num = 0;

  // First, try to increase till 2^14 = 16,384
  for (num = m_num_strides; num <= 16'384; num++) {
    if (len % num == 0)
      break;
  }

  if (len % num == 0) {
    m_num_strides = num;
    return;
  }

  // Second, try to decrease till 1, at which point it must work.
  for (num = m_num_strides; num > 0; num--) {
    if (len % num == 0)
      break;
  }

  m_num_strides = num;
}

auto sperr::Conditioner::condition(vecd_type& buf) -> std::pair<RTNType, meta_type>
{
  meta_type meta;
  meta.fill(0);
  double mean = 0.0;
  double rms = 0.0;

  // If divide by rms is requested, need to make sure rms is non-zero
  if (m_settings[1]) {
    if (std::all_of(buf.begin(), buf.end(), [](auto v) { return v == 0.0; }))
      return {RTNType::Error, meta};
  }

  m_adjust_strides(buf.size());

  //
  // The ordering of each operation is actually interesting. Between subtract
  // mean and divide by rms, it kind of makes sense to bring all values closer
  // to zero first and then do divide by rms, if you think about it...
  //

  // Perform subtract mean first
  if (m_settings[0]) {
    mean = m_calc_mean(buf);
    auto minus_mean = [mean](auto& v) { v -= mean; };
    std::for_each(buf.begin(), buf.end(), minus_mean);
    //std::cout<<mean<<std::endl;
  }
 
  // Perform divide_by_rms second
  if (m_settings[1]) {
    rms = m_calc_rms(buf);
    auto div_rms = [rms](auto& v) { v /= rms; };
    std::for_each(buf.begin(), buf.end(), div_rms);
  }

  // pack meta
  meta[0] = sperr::pack_8_booleans(m_settings);
  size_t pos = 1;
  std::memcpy(meta.data() + pos, &mean, sizeof(mean));
  pos += sizeof(mean);
  std::memcpy(meta.data() + pos, &rms, sizeof(rms));
  pos += sizeof(rms);  // NOLINT
  assert(pos == meta.size());

  return {RTNType::Good, meta};
}

auto sperr::Conditioner::inverse_condition(vecd_type& buf, const meta_type& meta) -> RTNType
{
  // unpack meta
  auto b8 = sperr::unpack_8_booleans(meta[0]);
  double mean = 0.0;
  double rms = 0.0;
  size_t pos = 1;
  std::memcpy(&mean, meta.data() + pos, sizeof(mean));
  pos += sizeof(mean);
  std::memcpy(&rms, meta.data() + pos, sizeof(rms));
  pos += sizeof(rms);  // NOLINT
  assert(pos == meta.size());

  m_adjust_strides(buf.size());

  // Perform inverse of divide_by_rms, which is multiply by rms
  if (b8[1]) {
    auto mul_rms = [rms](auto& v) { v *= rms; };
    std::for_each(buf.begin(), buf.end(), mul_rms);
  }

  // Perform inverse of subtract_mean, which is add mean.
  if (b8[0]) {
    auto plus_mean = [mean](auto& v) { v += mean; };
    std::for_each(buf.begin(), buf.end(), plus_mean);
  }

  return RTNType::Good;
}

auto sperr::Conditioner::test_constant(const vecd_type& buf) const -> std::pair<bool, meta_type>
{
  const uint64_t nval = buf.size();
  assert(nval > 0);

  const double val = buf[0];
  auto b8 = std::array<bool, 8>();
  b8.fill(false);

  if (std::all_of(buf.begin(), buf.end(), [val](auto v) { return v == val; }))
    b8[4] = true;

  // Prepare the meta block
  auto meta = meta_type();
  meta.fill(0);

  // It only makes sense to really prepare the meta block if the field tests to be constant.
  if (b8[4]) {
    // First byte of meta
    meta[0] = sperr::pack_8_booleans(b8);
    // Next 8 bytes of meta: the constant value
    size_t pos = 1;
    std::memcpy(meta.data() + pos, &val, sizeof(val));
    // Next 8 bytes of meta: how many constant values?
    pos += sizeof(val);
    std::memcpy(meta.data() + pos, &nval, sizeof(nval));
  }

  return {b8[4], meta};
}

auto sperr::Conditioner::parse_constant(const meta_type& meta) const
    -> std::tuple<bool, double, uint64_t>
{
  auto b8 = sperr::unpack_8_booleans(meta[0]);
  double val = 0.0;
  uint64_t nval = 0;

  // It only makes sense to really parse the meta block if the field tests to be constant.
  if (b8[4]) {
    // Next 8 bytes: the constant value
    size_t pos = 1;
    std::memcpy(&val, meta.data() + pos, sizeof(val));
    // Next 8 bytes: how many constant values:
    pos += sizeof(val);
    std::memcpy(&nval, meta.data() + pos, sizeof(nval));
  }

  return std::make_tuple(b8[4], val, nval);
}


#endif
