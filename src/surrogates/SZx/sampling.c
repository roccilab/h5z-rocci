#include "sampling.h"

void sample_1d_float(float* data, size_t blocksize, size_t sample_stride, size_t dims[], float** outBuffer, size_t* outSize, size_t sample_dims[]) {
    // Calculate the number of blocks
    size_t num_blocks = dims[0] / blocksize;

    // Calculate the size of the output buffer
    *outSize = ceil((float)num_blocks / sample_stride) * blocksize;
	size_t sampleBufferSize = *outSize;
	*outBuffer = malloc(sizeof(float) * (sampleBufferSize));

    // Initialize sample dimension
    sample_dims[0] = *outSize;
    size_t out_index = 0;

    // Iterate through blocks
    for (size_t i = 0; i < num_blocks; i+= sample_stride) {
        // Set the starting index of the current block
        size_t start_index = i * blocksize;

        // Iterate within the current block with the given stride
        for (size_t j = 0; j < blocksize; j ++) {

            // Calculate the index in the input data array
            size_t in_index = start_index + j;

            // Copy the data from the input to the output buffer
            (*outBuffer)[out_index] = data[in_index];
            out_index++;
        }
    }
}

void sample_2d_float(float* data, size_t blocksize, size_t sample_stride, size_t dims[], float** outBuffer, size_t* outSize, size_t sample_dims[]) {
    // Calculate the number of blocks in each dimension
    size_t blocksX = dims[0] / blocksize;
	size_t blocksY = dims[1] / blocksize;
    size_t num_blocks_x = ceil((float) blocksX / sample_stride);
    size_t num_blocks_y = ceil((float)blocksY / sample_stride);

    // Calculate the size of the output buffer
    *outSize = num_blocks_x * num_blocks_y * blocksize * blocksize;
	size_t sampleBufferSize = *outSize;
	*outBuffer = malloc(sizeof(float) * (sampleBufferSize));


    // Initialize sample dimensions
    sample_dims[0] = num_blocks_x * blocksize;
    sample_dims[1] = num_blocks_y * blocksize;


    // Iterate through blocks in each dimension with the given stride
    for (size_t i = 0; i < num_blocks_x; i++) {
        for (size_t j = 0; j < num_blocks_y; j++) {
            // Set the starting index of the current block
            size_t start_index_x = (i * sample_stride) * blocksize;
            size_t start_index_y = (j * sample_stride) * blocksize;

            // Iterate within the current block with the given stride
            for (size_t x = 0; x < blocksize; x++) {
                for (size_t y = 0; y < blocksize; y++) {
                    // Calculate the index in the input data array
                    size_t in_index_x = start_index_x + x;
                    size_t in_index_y = start_index_y + y;
                    size_t data_index = in_index_x * dims[1] + in_index_y;

                    // Calculate the output index using a linear combination
                    size_t out_index = ((i * blocksize + x) * sample_dims[1]) + (j * blocksize + y);

                    // Copy the data from the input to the output buffer
                    (*outBuffer)[out_index] = data[data_index];
                }
            }
        }
    }
}


void sample_3d_float(float* data, size_t blocksize, size_t sample_stride, size_t dims[], float** outBuffer, size_t* outSize, size_t sample_dims[]) {
    // Calculate the number of blocks in each dimension
    size_t blocksX = dims[0] / blocksize;
    size_t blocksY = dims[1] / blocksize;
    size_t blocksZ = dims[2] / blocksize;
    size_t num_blocks_x = ceil((float) blocksX / sample_stride);
    size_t num_blocks_y = ceil((float) blocksY / sample_stride);
    size_t num_blocks_z = ceil((float) blocksZ / sample_stride);

    // Calculate the size of the output buffer
    *outSize = num_blocks_x * num_blocks_y * num_blocks_z * blocksize * blocksize * blocksize;
	size_t sampleBufferSize = *outSize;
	*outBuffer = malloc(sizeof(float) * (sampleBufferSize));

    // Initialize sample dimensions
    sample_dims[0] = num_blocks_x * blocksize;
    sample_dims[1] = num_blocks_y * blocksize;
    sample_dims[2] = num_blocks_z * blocksize;

    // Iterate through blocks in each dimension with the given stride
    for (size_t i = 0; i < num_blocks_x; i++) {
        for (size_t j = 0; j < num_blocks_y; j++) {
            for (size_t k = 0; k < num_blocks_z; k++) {
                // Set the starting index of the current block
                size_t start_index_x = (i * sample_stride) * blocksize;
                size_t start_index_y = (j * sample_stride) * blocksize;
                size_t start_index_z = (k * sample_stride) * blocksize;

                // Iterate within the current block with the given stride
                for (size_t x = 0; x < blocksize; x++) {
                    for (size_t y = 0; y < blocksize; y++) {
                        for (size_t z = 0; z < blocksize; z++) {
                            // Calculate the index in the input data array
                            size_t in_index_x = start_index_x + x;
                            size_t in_index_y = start_index_y + y;
                            size_t in_index_z = start_index_z + z;
                            size_t data_index = (in_index_x * dims[1] * dims[2]) + (in_index_y * dims[2]) + in_index_z;

                            // Calculate the output index using a linear combination
                            size_t out_index = (((i * blocksize + x) * sample_dims[1] * sample_dims[2]) +
                                                ((j * blocksize + y) * sample_dims[2]) +
                                                (k * blocksize + z));

                            // Copy the data from the input to the output buffer
                            (*outBuffer)[out_index] = data[data_index];
                        }
                    }
                }
            }
        }
    }
}