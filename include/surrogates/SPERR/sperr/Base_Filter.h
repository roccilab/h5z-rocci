#ifndef BASE_FILTER_H
#define BASE_FILTER_H

//
// This is the base class for all custom filters, which does not perform any meaningful task.
//

#include "sperr_helper.h"

namespace sperr {

class Base_Filter {
 public:
  virtual ~Base_Filter() = default;

  //
  // Action items
  //
  // Given data `buf` and its dimensions `dims`, apply the filter.
  // The filter can create a header that facilitates an inverse filtering step
  // and return that header.
  // Note that the header should always save its length.
  virtual auto apply_filter(vecd_type& buf, dims_type dims) -> vec8_type;

  // Given data `buf`, its dimensions `dims`, and a header, apply the inverse filter.
  // Return true if the operation is successful; false otherwise.
  virtual auto inverse_filter(vecd_type& buf, dims_type dims, const void* header, size_t header_len)
      -> bool;

  // Given a header, retrieve its length.
  virtual auto header_size(const void* header) const -> size_t;
};

};  // namespace sperr
auto sperr::Base_Filter::apply_filter(vecd_type& buf, dims_type dims) -> vec8_type
{
  auto empty = vec8_type();
  return empty;
}

auto sperr::Base_Filter::inverse_filter(vecd_type& buf,
                                        dims_type dims,
                                        const void* header,
                                        size_t header_len) -> bool
{
  return true;
}

auto sperr::Base_Filter::header_size(const void* header) const -> size_t
{
  return 0;
}
#endif