#ifndef SPERR3D_COMPRESSOR_H
#define SPERR3D_COMPRESSOR_H

#include "CDF97.h"
#include "Conditioner.h"
#include "SPECK3D.h"
#include "SPERR.h"

//#ifdef USE_ZSTD
#include "zstd.h"
//#endif

#include <algorithm>
#include <cassert>
#include <cstring>

#include<iostream>
namespace sperr {

class SPERR3D_Compressor {
 public:
  // Accept incoming data: copy from a raw memory block
  template <typename T>
  auto copy_data(const T* p,      // Input: pointer to the memory block
                 size_t len,      // Input: number of values
                 dims_type dims)  // Input: dimension of the 3D volume
      -> RTNType;

  // Accept incoming data: take ownership of a memory block
  auto take_data(std::vector<double>&& buf, dims_type dims) -> RTNType;

  //void toggle_conditioning(Conditioner::settings_type);

  auto set_comp_params(size_t bit_budget, double psnr, double pwe) -> RTNType;

  // Return 1) the number of outliers, and 2) the number of bytes to encode them.
  auto get_outlier_stats() const -> std::pair<size_t, size_t>;

  auto compress() -> RTNType;

  auto view_encoded_bitstream() const -> const std::vector<uint8_t>&;
  auto release_encoded_bitstream() -> std::vector<uint8_t>&&;
  void set_eb_coeff(const double & coeff);

  void set_skip_wave(const bool & skip);

 private:
  dims_type m_dims = {0, 0, 0};
  vecd_type m_val_buf;

  Conditioner m_conditioner;
  CDF97 m_cdf;
  SPECK3D m_encoder;

 // Conditioner::settings_type m_conditioning_settings = {true, false, false, false};

  // Store bitstreams from the conditioner and SPECK encoding, and the overall bitstream.
 // Conditioner::meta_type m_condi_stream;
  vec8_type m_condi_stream;
  vec8_type m_encoded_stream;

  size_t m_bit_budget = 0;  // Total bit budget, including headers etc.
  double m_target_psnr = sperr::max_d;
  double m_target_pwe = 0.0;
  double eb_coeff=1.5;
  bool skip_wave=false;

  SPERR m_sperr;
  vec8_type m_sperr_stream;
  vecd_type m_val_buf2;        // Copy of `m_val_buf` that's used for outlier coding.
  std::vector<Outlier> m_LOS;  // List of OutlierS

//#ifdef USE_ZSTD
  vec8_type m_zstd_buf;
  std::unique_ptr<ZSTD_CCtx, decltype(&ZSTD_freeCCtx)> m_cctx = {nullptr, &ZSTD_freeCCtx};
//#endif

  auto m_assemble_encoded_bitstream() -> RTNType;
};

};  // namespace sperr


void sperr::SPERR3D_Compressor::set_skip_wave(const bool & skip){
  skip_wave=skip;
}


void sperr::SPERR3D_Compressor::set_eb_coeff(const double & coeff){
  eb_coeff=coeff;
}

template <typename T>
auto sperr::SPERR3D_Compressor::copy_data(const T* p, size_t len, sperr::dims_type dims) -> RTNType
{
  static_assert(std::is_floating_point<T>::value, "!! Only floating point values are supported !!");

  if (len != dims[0] * dims[1] * dims[2])
    return RTNType::WrongDims;

  m_val_buf.resize(len);
  std::copy(p, p + len, m_val_buf.begin());

  m_dims = dims;

  return RTNType::Good;
}
template auto sperr::SPERR3D_Compressor::copy_data(const double*, size_t, dims_type) -> RTNType;
template auto sperr::SPERR3D_Compressor::copy_data(const float*, size_t, dims_type) -> RTNType;

auto sperr::SPERR3D_Compressor::take_data(sperr::vecd_type&& buf, sperr::dims_type dims) -> RTNType
{
  if (buf.size() != dims[0] * dims[1] * dims[2])
    return RTNType::WrongDims;

  m_val_buf = std::move(buf);
  m_dims = dims;

  return RTNType::Good;
}

auto sperr::SPERR3D_Compressor::view_encoded_bitstream() const -> const std::vector<uint8_t>&
{
  return m_encoded_stream;
}

auto sperr::SPERR3D_Compressor::release_encoded_bitstream() -> std::vector<uint8_t>&&
{
  return std::move(m_encoded_stream);
}

auto sperr::SPERR3D_Compressor::compress() -> RTNType
{
  const auto total_vals = m_dims[0] * m_dims[1] * m_dims[2];
  if (m_val_buf.empty() || m_val_buf.size() != total_vals)
    return RTNType::Error;
  //m_condi_stream.fill(0);
  m_condi_stream.clear();
  m_encoded_stream.clear();

  m_sperr_stream.clear();
  m_val_buf2.clear();
  m_LOS.clear();
  //std::cout<<"p1"<<std::endl;
  /*
  // Believe it or not, there are constant fields passed in for compression!
  // Let's detect that case and skip the rest of the compression routine if it occurs.
  auto constant = m_conditioner.test_constant(m_val_buf);
  if (constant.first) {
    m_condi_stream = constant.second;
    auto tmp = m_assemble_encoded_bitstream();
    return tmp;
  }
  */
  // Keep track of data range before and after the conditioning step, in case they change.
  // This is only used in `FixedPWE` mode though.
  auto range_before = double{0.0};
  auto range_after = double{0.0};


  // Find out the compression mode, and initialize data members accordingly.
  const auto mode = sperr::compression_mode(m_bit_budget, m_target_psnr, m_target_pwe);
  assert(mode != CompMode::Unknown);
  /*
  if (mode == sperr::CompMode::FixedPSNR) {
    // Calculate the original data range and pass it to the encoder.
    auto [min, max] = std::minmax_element(m_val_buf.cbegin(), m_val_buf.cend());
    auto range = *max - *min;
    m_encoder.set_data_range(range);
  }
  else if (mode == sperr::CompMode::FixedPWE) {
    */ 
  if (mode == sperr::CompMode::FixedPWE ) {
    // Make a copy of the original data for outlier correction use.
    m_val_buf2.resize(total_vals);
    std::copy(m_val_buf.cbegin(), m_val_buf.cend(), m_val_buf2.begin());
    // auto [min, max] = std::minmax_element(m_val_buf.cbegin(), m_val_buf.cend());//commented for saving time as currently custom filter is not used.
    //range_before = *max - *min;
  }
  //std::cout<<"p2"<<std::endl;
  // Step 1: data goes through the conditioner
  if(1){
    /*
    m_conditioner.toggle_all_settings(m_conditioning_settings);
    auto [rtn, condi_meta] = m_conditioner.condition(m_val_buf);
    if (rtn != RTNType::Good)
    */
    m_condi_stream = m_conditioner.condition(m_val_buf, m_dims);
    assert(!m_condi_stream.empty());
    // Step 1.1: believe it or not, there are constant fields passed in for compression!
    // Let's detect that case and skip the rest of the compression routine if it occurs.
    if (m_conditioner.is_constant(m_condi_stream[0])) {
      auto rtn = m_assemble_encoded_bitstream();
      return rtn;
    //m_condi_stream = condi_meta;
    }
    if (mode == sperr::CompMode::FixedPSNR) {
      // Calculate data range using the conditioned data, and pass it to the encoder.
      auto [min, max] = std::minmax_element(m_val_buf.cbegin(), m_val_buf.cend());
      auto range = *max - *min;
      m_encoder.set_data_range(range);
    }
    else if (mode == sperr::CompMode::FixedPWE &&
             m_conditioner.has_custom_filter(m_condi_stream[0])) {
      // Only re-calculate data range when there's custom filter enabled in the conditioner.
      auto [min, max] = std::minmax_element(m_val_buf.cbegin(), m_val_buf.cend());
      range_after = *max - *min;
    }

    //std::cout<<"p3"<<std::endl;


    // Step 2: wavelet transform
    //rtn = m_cdf.take_data(std::move(m_val_buf), m_dims);
    auto rtn = m_cdf.take_data(std::move(m_val_buf), m_dims);
    if (rtn != RTNType::Good)
      return rtn;
    m_cdf.dwt3d();

    // Step 3: SPECK encoding
    const auto & coeffs=m_cdf.release_data();

    //sperr::write_n_bytes("sperr.dwt",coeffs.size()*sizeof(double),coeffs.data());

    rtn = m_encoder.take_data(m_cdf.release_data(), m_dims);
    if (rtn != RTNType::Good)
      return rtn;

    //std::cout<<"p4"<<std::endl;

  }
  else{

    // Step 1: Direct SPECK encoding
    //const auto & coeffs=m_cdf.release_data();

    //sperr::write_n_bytes("sperr.dwt",coeffs.size()*sizeof(double),coeffs.data());
    const auto header_size = 1 + sizeof(double) + 0;//0 is the custom filter header size.
    m_condi_stream.resize(header_size,0);
    auto b8=sperr::unpack_8_booleans(m_condi_stream[0]);
    b8[2]=true;
    m_condi_stream[0]=sperr::pack_8_booleans(b8);
    auto rtn = m_encoder.take_data(std::move(m_val_buf), m_dims);
    if (rtn != RTNType::Good)
      return rtn;

  }

  

  auto speck_budget = size_t{0};
  if (m_bit_budget == sperr::max_size)
    speck_budget = sperr::max_size;
  else
    speck_budget = m_bit_budget - m_condi_stream.size() * 8;
  m_encoder.set_eb_coeff(eb_coeff);
  //auto rtn = m_encoder.set_comp_params(speck_budget, m_target_psnr, m_target_pwe);
  // In the FixedPWE mode, in case there's custom filter, we scale the PWE tolerance
  auto speck_pwe = m_target_pwe;
  if (mode == sperr::CompMode::FixedPWE && m_conditioner.has_custom_filter(m_condi_stream[0])) {
    assert(range_before != 0.0);
    speck_pwe *= range_after / range_before;
  }

  auto rtn = m_encoder.set_comp_params(speck_budget, m_target_psnr, speck_pwe);
  if (rtn != RTNType::Good)
    return rtn;

  rtn = m_encoder.encode();
 //  std::cout<<"p5"<<std::endl;
  if (rtn != RTNType::Good)
    return rtn;
  //std::cout<<m_encoder.count_LSP()<<std::endl;
  // We can copy the encoded speck stream to `m_encoded_stream` at the appropriate position.
  const auto& speck_stream = m_encoder.view_encoded_bitstream();
 // if (speck_stream.empty())
  //  return RTNType::Error;
  assert(!speck_stream.empty());
  const size_t condi_speck_len = m_condi_stream.size() + speck_stream.size();
  m_encoded_stream.resize(condi_speck_len);
  std::copy(speck_stream.cbegin(), speck_stream.cend(),
            m_encoded_stream.begin() + m_condi_stream.size());

  // Step 4: Outlier correction if in FixedPWE mode and not using speck only.
  if (mode == sperr::CompMode::FixedPWE and !skip_wave) {
    // Step 4.1: IDWT using quantized coefficients to have a reconstruction.
    auto qz_coeff = m_encoder.release_quantized_coeff();

    //sperr::write_n_bytes("sperr.dwt.cmpdec",qz_coeff.size()*sizeof(double),qz_coeff.data());


    assert(!qz_coeff.empty());
    m_cdf.take_data(std::move(qz_coeff), m_dims);
    m_cdf.idwt3d();
    m_val_buf = m_cdf.release_data();
    //m_conditioner.inverse_condition(m_val_buf, m_condi_stream);
    m_conditioner.inverse_condition(m_val_buf, m_dims, m_condi_stream);
   //  std::cout<<"p5.1"<<std::endl;
    // Step 4.2: Find all outliers
    for (size_t i = 0; i < total_vals; i++) {
      const auto diff = m_val_buf2[i] - m_val_buf[i];
      if (std::abs(diff) >= m_target_pwe)
        m_LOS.emplace_back(i, diff);
    }

    // Step 4.3: Code located outliers
    if (!m_LOS.empty()) {
      m_sperr.set_tolerance(m_target_pwe);
      m_sperr.set_length(total_vals);
      m_sperr.copy_outlier_list(m_LOS);
      rtn = m_sperr.encode();
      if (rtn != RTNType::Good)
        return rtn;
      m_sperr_stream = m_sperr.get_encoded_bitstream();
      if (m_sperr_stream.empty())
        return RTNType::Error;
    }
  }
  // std::cout<<"p6"<<std::endl;
  rtn = m_assemble_encoded_bitstream();

  return rtn;
}

auto sperr::SPERR3D_Compressor::m_assemble_encoded_bitstream() -> RTNType
{
  // Copy over the condi stream.
  // `m_encoded_stream` is either empty, or is already allocated space and filled with the speck
  // stream.
  if (m_encoded_stream.empty())
    m_encoded_stream.resize(m_condi_stream.size());
  std::copy(m_condi_stream.begin(), m_condi_stream.end(), m_encoded_stream.begin());

  // If there's outlier correction stream, copy it over too.
  if (!m_sperr_stream.empty()) {
    const auto condi_speck_len = m_encoded_stream.size();
    assert(condi_speck_len > m_condi_stream.size());
    m_encoded_stream.resize(condi_speck_len + m_sperr_stream.size());
    std::copy(m_sperr_stream.cbegin(), m_sperr_stream.cend(),
              m_encoded_stream.begin() + condi_speck_len);
  }

//#ifdef USE_ZSTD
  if (m_cctx == nullptr) {
    auto* ctx_p = ZSTD_createCCtx();
    if (ctx_p == nullptr)
      return RTNType::ZSTDError;
    else
      m_cctx.reset(ctx_p);
  }

  const size_t comp_upper_size = ZSTD_compressBound(m_encoded_stream.size());
  m_zstd_buf.resize(comp_upper_size);
  const size_t comp_size =
      ZSTD_compressCCtx(m_cctx.get(), m_zstd_buf.data(), comp_upper_size, m_encoded_stream.data(),
                        m_encoded_stream.size(), ZSTD_CLEVEL_DEFAULT + 6);
  if (ZSTD_isError(comp_size))
    return RTNType::ZSTDError;
  else {
    // Note: when the encoded stream is only a few kilobytes or smaller, the ZSTD compressed
    //       output can be larger.
    m_encoded_stream.resize(comp_size);
    std::copy(m_zstd_buf.cbegin(), m_zstd_buf.cbegin() + comp_size, m_encoded_stream.begin());
  }
//#endif

  return RTNType::Good;
}

auto sperr::SPERR3D_Compressor::get_outlier_stats() const -> std::pair<size_t, size_t>
{
  return {m_LOS.size(), m_sperr_stream.size()};
}

auto sperr::SPERR3D_Compressor::set_comp_params(size_t budget, double psnr, double pwe) -> RTNType
{
  // First set those ones that only need a plain copy
  m_target_psnr = psnr;
  m_target_pwe = pwe;

  // Second set bit budget, which would require a little manipulation.
  if (budget == sperr::max_size) {
    m_bit_budget = sperr::max_size;
    return RTNType::Good;
  }

  if (budget < m_condi_stream.size() * 8)
    return RTNType::InvalidParam;
  else {
    m_bit_budget = budget;
    return RTNType::Good;
  }
}
/*
void sperr::SPERR3D_Compressor::toggle_conditioning(sperr::Conditioner::settings_type b4)
{
  m_conditioning_settings = b4;
}

*/

#endif
