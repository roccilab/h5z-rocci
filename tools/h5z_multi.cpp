#include <fstream>
#include <iterator>
#include <memory>
#include <vector>
#include <cmath>

#include "hdf5.h"
#include "H5PLextern.h"
#include "debug.hpp"
#include "rocci_config.hpp"
#include "h5z_multi.hpp"
#include "select.hpp"

hid_t H5Z_ROCCI_ERRCLASS = -1;

__attribute__((constructor))
void filter_debug_constructor() {
    fprintf(stderr, "Filter plugin loaded successfully!\n");
}

// filter definition
__attribute__((used))
__attribute__((visibility("default")))
const H5Z_class2_t H5Z_ROCCI[1] = {{
    H5Z_CLASS_T_VERS,                                       /* H5Z_class_t version */
    H5Z_FILTER_ROCCI,                                       /* Filter id number */
    1,                                                      /* encoder_present flag (set to true) */
    1,                                                      /* decoder_present flag (set to true) */
    "TEST ROCCI",                                           /* Filter name for debugging */
    NULL,                                                   /* The "can apply" callback */
    H5Z_ROCCI_set_local,                                    /* The "set local" callback */
    static_cast<H5Z_func_t>(H5Z_filter_ROCCI),              /* The actual filter function */
}};

H5PL_type_t H5PLget_plugin_type(void) { return H5PL_TYPE_FILTER; }

const void *H5PLget_plugin_info(void) { return H5Z_ROCCI; }

herr_t set_ROCCI_conf_to_H5(const hid_t propertyList, Config &conf) {
    static char const *_funcname_ = "set_ROCCI_conf_to_H5";

    // save conf into cd_values
    size_t cd_nelmts = std::ceil(conf.size_est() / 1.0 / sizeof(int));
    std::vector<unsigned int> cd_values(cd_nelmts, 0);
    auto buffer = reinterpret_cast<unsigned char *>(cd_values.data());

    auto confSizeReal = conf.save(buffer);
    cd_nelmts = std::ceil(confSizeReal / 1.0 / sizeof(int));

    auto filter = H5Pget_filter_by_id(propertyList, H5Z_FILTER_ROCCI, H5Z_FLAG_MANDATORY, NULL, NULL, 0, NULL,
                                        NULL);  // check if filter is set
    if (0 > filter) {  // filter not set, set filter. Notice that calling H5Pset_filter twice with the same filter id
                         // will cause unexpected errors for decompression
        if (0 > H5Pset_filter(propertyList, H5Z_FILTER_ROCCI, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values.data())) {
            H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
        }
    } else {  // filter already set, update filter
        if (0 > H5Pmodify_filter(propertyList, H5Z_FILTER_ROCCI, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values.data()))
            H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
    }

    return 1;
}

herr_t get_ROCCI_conf_from_H5(const hid_t propertyList, Config &conf) {
    static char const *_funcname_ = "get_ROCCI_conf_from_H5";

    size_t cd_nelmts = std::ceil(conf.size_est() / 1.0 / sizeof(int));
    std::vector<unsigned int> cd_values(cd_nelmts, 0);

    // read cd_values from HDF5
    // note that cd_nelmts must be non-zero, otherwise, cd_values cannot be filled.
    if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_ROCCI, H5Z_FLAG_MANDATORY, &cd_nelmts, cd_values.data(), 0, NULL,
                                NULL))
        H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current ROCCI cd_values");

    // load cd_values into config
    if (cd_nelmts != 0) {
        auto buffer = reinterpret_cast<const unsigned char *>(cd_values.data());
        conf.load(buffer);
    }
    return 1;
}

static herr_t H5Z_ROCCI_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {
    std::cout << "start H5Z_ROCCI_set_local" << std::endl;

    static char const *_funcname_ = "H5Z_ROCCI_set_local";

    Config conf;
    get_ROCCI_conf_from_H5(dcpl_id, conf);
    
    // read datatype and dims from HDF5
    H5T_class_t dclass;
    if (0 > (dclass = H5Tget_class(type_id))) H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");

    size_t dsize;
    if (0 == (dsize = H5Tget_size(type_id))) H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");

    int ndims;
    hsize_t dims_all[H5S_MAX_RANK];
    if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims_all, NULL)))
        H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");
    std::vector<size_t> dims(dims_all, dims_all + ndims);

    if (dclass == H5T_FLOAT)
        conf.dataType = dsize == 4 ? ROCCI_FLOAT : ROCCI_DOUBLE;
    else if (dclass == H5T_INTEGER) {
        H5T_sign_t dsign;
        if (0 > (dsign = H5Tget_sign(type_id)))
            H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");
        if (dsign == H5T_SGN_NONE)  // unsigned
        {
            switch (dsize) {
                case 1:
                    conf.dataType = ROCCI_UINT8;
                    break;
                case 2:
                    conf.dataType = ROCCI_UINT16;
                    break;
                case 4:
                    conf.dataType = ROCCI_UINT32;
                    break;
                case 8:
                    conf.dataType = ROCCI_UINT64;
                    break;
            }
        } else {
            switch (dsize) {
                case 1:
                    conf.dataType = ROCCI_INT8;
                    break;
                case 2:
                    conf.dataType = ROCCI_INT16;
                    break;
                case 4:
                    conf.dataType = ROCCI_INT32;
                    break;
                case 8:
                    conf.dataType = ROCCI_INT64;
                    break;
            }
        }
    } else {
        H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
    }

    conf.setDims(dims.begin(), dims.end());
    set_ROCCI_conf_to_H5(dcpl_id, conf);

    return 1;
}

template <typename T>
size_t ROCCI_compress(Config &conf, void *data_void_ptr, char *&compressedData) {
    // T *data = static_cast<T *>(data_void_ptr);
    // size_t nbytes = conf.num * sizeof(T);
    // compressedData = static_cast<char *>(malloc(nbytes));
    // memcpy(compressedData, data, nbytes);
    // return nbytes;

    size_t compressed_size;
    compressed_size = ROCCI_select_compress(conf, static_cast<T *>(data_void_ptr), compressedData, compressed_size);
    return compressed_size;
}

template <typename T>
void ROCCI_decompress(Config &conf, char *compressedData, void *decompressedData) {
    // Placeholder for decompression logic
    // This function should decompress the data based on the configuration in conf
    // and store the result in decompressedData.
    // For now, we just simulate decompression by copying data.

    ROCCI_select_decompress(conf, static_cast<T *>(decompressedData), compressedData);
}

template <typename T>
void process_data(Config &conf, void **buf, size_t *buf_size, size_t nbytes, bool is_decompress) {
    if (is_decompress) {
        T *processedData = static_cast<T *>(malloc(conf.num * sizeof(T)));
        ROCCI_decompress<T>(conf, static_cast<char *>(*buf), processedData);
        free(*buf);
        *buf = processedData;
        *buf_size = conf.num * sizeof(T);
    } else {
        char *cmpData = nullptr;
        *buf_size = ROCCI_compress<T>(conf, static_cast<T *>(*buf), cmpData);
        free(*buf);
        *buf = cmpData;
    }
}

/**
 * https://docs.hdfgroup.org/hdf5/v1_14/_f_i_l_t_e_r.html
 * The flags, cd_nelmts, and cd_values are the same as for the H5Pset_filter() function with the additional flag
 * H5Z_FLAG_REVERSE which is set when the filter is called as part of the input pipeline. The input buffer is pointed to
 * by *buf and has a total size of *buf_size bytes but only nbytes are valid data. The filter should perform the
 * transformation in place if possible and return the number of valid bytes or zero for failure. If the transformation
 * cannot be done in place then the filter should allocate a new buffer with malloc() and assign it to *buf, assigning
 * the allocated size of that buffer to *buf_size. The old buffer should be freed by calling free().
 */
static size_t H5Z_filter_ROCCI(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes,
                                size_t *buf_size, void **buf) {
    std::cout << "start H5Z_filter_ROCCI" << std::endl;

    if (cd_nelmts == 0)  // this is special data such as string, which should not be treated as values.
        return nbytes;

    Config conf;
    auto buffer = reinterpret_cast<const unsigned char *>(cd_values);
    conf.load(buffer);

    bool is_decompress = flags & H5Z_FLAG_REVERSE;
    // get data type and dimensions from set_local

    // if (is_decompress) printf("doing decompression\n");
    // else printf("doing compression\n");
    // fflush(stdout);
    
    switch (conf.dataType) {
        case ROCCI_FLOAT:
            process_data<float>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_DOUBLE:
            process_data<double>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_INT8:
            process_data<int8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_UINT8:
            process_data<uint8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_INT16:
            process_data<int16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_UINT16:
            process_data<uint16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_INT32:
            process_data<int32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_UINT32:
            process_data<uint32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_INT64:
            process_data<int64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case ROCCI_UINT64:
            process_data<uint64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        default:
            std::cerr << (is_decompress ? "Decompression" : "Compression") << " Error: Unknown Datatype" << std::endl;
            std::exit(EXIT_FAILURE);
    }

    return *buf_size;
}
