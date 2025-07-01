#ifndef __ROCCI_CONFIG_HPP__
#define __ROCCI_CONFIG_HPP__

#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <cstring>
#include <algorithm>

#define ROCCI_FLOAT 0
#define ROCCI_DOUBLE 1
#define ROCCI_UINT8 2
#define ROCCI_INT8 3
#define ROCCI_UINT16 4
#define ROCCI_INT16 5
#define ROCCI_UINT32 6
#define ROCCI_INT32 7
#define ROCCI_UINT64 8
#define ROCCI_INT64 9

enum ALGO { ALGO_SZ3, ALGO_ZFP, ALGO_SZx, ALGO_SPERR };
constexpr const char *ALGO_STR[] = { "SZ3", "ZFP", "SZx", "SPERR" };

constexpr const ALGO ALGO_OPTIONS[] = {  
    ALGO::ALGO_SZ3,
    ALGO::ALGO_ZFP,
    ALGO::ALGO_SZx,
    ALGO::ALGO_SPERR
};

template <class T>
const char *enum2Str(T e) {
    if (std::is_same<T, ALGO>::value) {
        return ALGO_STR[e];
    } else {
        fprintf(stderr, "invalid enum type for enum2Str()\n ");
        throw std::invalid_argument("invalid enum type for enum2Str()");
    }
}

class Config {
   public:
    template <class... Dims>
    Config(Dims... args) {
        dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
        setDims(dims.begin(), dims.end());
    }

    template <class Iter>
    size_t setDims(Iter begin, Iter end) {
        auto dims_ = std::vector<size_t>(begin, end);
        dims.clear();
        for (auto dim : dims_) {
            if (dim > 1) {
                dims.push_back(dim);
            }
        }
        if (dims.empty()) {
            dims = {1};
        }
        N = dims.size();
        num = std::accumulate(dims.begin(), dims.end(), static_cast<size_t>(1), std::multiplies<size_t>());
        return num;
    }

    size_t save(unsigned char *&c) const {
        auto c0 = c;
        write(N, c);
        // write(dims.data(), dims.size(), c);
        auto bitWidth = vector_bit_width(dims);
        write(bitWidth, c);
        vector2bytes(dims, bitWidth, c);

        write(num, c);
        write(cmprAlgo, c);

        write(absErrorBound, c);

        write(dataType, c);

        return c - c0;
    }

    void load(const unsigned char *&c) {
        read(N, c);

        uint8_t bitWidth;
        read(bitWidth, c);
        dims = bytes2vector<size_t>(c, bitWidth, N);
        read(num, c);
        read(cmprAlgo, c);
        
        read(absErrorBound, c);

        read(dataType, c);
    }

    size_t size_est() {
        return sizeof(Config) + sizeof(size_t) * 32; // 32 for the dims
    }

    char N = 0;
    std::vector<size_t> dims;
    size_t num = 0;
    uint8_t cmprAlgo = ALGO_SZ3;
    double absErrorBound = 1e-1;
    uint8_t dataType = ROCCI_FLOAT;  // dataType is only used in HDF5 filter

    private:

    using uchar = unsigned char;

    template <class T1>
    void write(T1 const var, uchar *&compressed_data_pos) const {
        memcpy(compressed_data_pos, &var, sizeof(T1));
        compressed_data_pos += sizeof(T1);
    }

    // read variable
    template <class T1>
    void read(T1 &var, uchar const *&compressed_data_pos) {
        memcpy(&var, compressed_data_pos, sizeof(T1));
        compressed_data_pos += sizeof(T1);
    }
    
    template <typename T>
    uint8_t vector_bit_width(const std::vector<T> &data) const {
        if (data.empty()) return 0;
        T max_value = *std::max_element(data.begin(), data.end());
        uint8_t bits = 0;
        while (max_value > 0) {
            max_value >>= 1;
            ++bits;
        }
        return bits;
    }

    template <typename T>
    void vector2bytes(const std::vector<T> &data, uint8_t bit_width, unsigned char *&c) const {
        if (data.empty()) return;

        size_t current_bit = 0;
        size_t byte_index = 0;
        unsigned char current_byte = 0;

        for (T value : data) {
            size_t bits_remaining = bit_width;
            while (bits_remaining > 0) {
                size_t space_in_current_byte = 8 - (current_bit % 8);
                size_t bits_to_write = std::min(bits_remaining, space_in_current_byte);
                size_t bits_shift = (bit_width - bits_remaining);
                unsigned char bits_to_store = (value >> bits_shift) & ((1 << bits_to_write) - 1);

                current_byte |= (bits_to_store << (current_bit % 8));
                current_bit += bits_to_write;
                bits_remaining -= bits_to_write;

                if (current_bit % 8 == 0) {
                    c[byte_index++] = current_byte;
                    current_byte = 0;
                }
            }
        }

        if (current_bit % 8 != 0) {
            c[byte_index++] = current_byte;
        }

        c += byte_index;
    }

    template <typename T>
    std::vector<T> bytes2vector(const unsigned char *&c, uint8_t bit_width, size_t num_elements) {
        // uint8_t bit_width = *c++;

        std::vector<T> data(num_elements);

        size_t total_bits = num_elements * bit_width;
        size_t total_bytes = (total_bits + 7) / 8;

        for (size_t i = 0; i < num_elements; ++i) {
            T value = 0;
            for (uint8_t j = 0; j < bit_width; ++j) {
                size_t bit_index = i * bit_width + j;
                size_t byte_index = bit_index / 8;
                size_t bit_offset = bit_index % 8;

                value |= ((c[byte_index] >> bit_offset) & 1) << j;
            }
            data[i] = value;
        }

        c += total_bytes;

        return data;
    }
};

#endif
