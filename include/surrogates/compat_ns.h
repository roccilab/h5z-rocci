#pragma once
// Compatibility shim for libpressio 1.x namespace changes
#include <libpressio_ext/cpp/compressor.h>
namespace libpressio {
  using compressors::libpressio_compressor_plugin;
}
