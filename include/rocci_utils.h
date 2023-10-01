#include <stdlib.h>
#include <rocci.h>

size_t computeNbEle(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

enum pressio_dtype rocci_dtype_to_pressio_dtype(int dataType);

struct pressio_compressor* get_compressor(ROCCI_Setting rocci_config, struct pressio* library);

struct pressio_options* get_compressor_options(ROCCI_Setting rocci_config);

struct pressio_data* create_pressio_data(int dataType, void* data, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);