#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <SZp/SZp_float.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#ifdef _OPENMP
#include "omp.h"
#endif

static inline int quantization_f32(float data, float recipPrecision)
{
    float dataRecip = data*recipPrecision;
    int s = dataRecip>=-0.5f?0:1;
    return (int)(dataRecip+0.5f) - s;
}

// compute size of sign byte array
size_t estimate_signbytelength(size_t intArrayLength) {
    size_t byteLength = 0;
    if (intArrayLength % 8 == 0)
		byteLength = intArrayLength / 8;
	else
		byteLength = intArrayLength / 8 + 1;
    return byteLength;
}

// compute size of saved bits array
size_t estimate_savedbitsbytelength(size_t intArrayLength, unsigned int bit_count)
{
    unsigned int byte_count = 0;
	unsigned int remainder_bit = 0;
	byte_count = bit_count / 8;
	remainder_bit = bit_count % 8;
	size_t byteLength = byte_count * intArrayLength + (remainder_bit * intArrayLength - 1) / 8 + 1;
	size_t byte_offset = byte_count * intArrayLength;
	if (remainder_bit == 0)
	{
		byteLength = byte_offset;
	}
    return byteLength;
}

double SZp_float_estimate_compress(float* data, size_t nbEle, int blockSize, float absErrBound, size_t sample_stride)
{
    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int new_block_size = block_size - 1;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;

    size_t nbSampleEle = 0;

    unsigned int outSize = 0;

    #pragma omp parallel
    {
        #pragma omp single
        {
            nbThreads = omp_get_num_threads();
            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
            remainder = nbEle % nbThreads;
            num_full_block_in_tb = (threadblocksize) / block_size; // minus 1 because the first one will be direcrly saved
            num_remainder_in_tb = (threadblocksize) % block_size;
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
        }

        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        size_t outSize_perthread = 0;
        // qutization centered with 0
        int tid = omp_get_thread_num();
        int lo = tid * threadblocksize;
        int hi = (tid + 1) * threadblocksize;

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned int signbytelength = 0; // save the bytelength
        unsigned int savedbitsbytelength = 0;

        if (num_full_block_in_tb > 0)
        {
            // iterate blocks by sample stride
            for (i = lo; i < hi - num_remainder_in_tb; i = i + (block_size * sample_stride))
            {
                max = 0;
                prior = (data[i]) * inver_bound;
                outSize_perthread += sizeof(unsigned int);

                for (j = 0; j < new_block_size; j++)
                {
                    current = (data[i + j + 1]) * inver_bound;
                    nbSampleEle++;

                    diff = current - prior;
                    prior = current;
                    if (diff > 0)
                    {
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                    else if (diff < 0)
                    {
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                }
                if (max == 0) // can be memset
                {
                    // only count the fixed length
                    outSize_perthread++;
                }
                else
                {
                    // count the fixed length, sign bits, and the saved bits
                    bit_count = (int)(log2f(max)) + 1;
                    outSize_perthread++;

                    signbytelength = estimate_signbytelength(new_block_size);
                    outSize_perthread += signbytelength;

                    savedbitsbytelength = estimate_savedbitsbytelength(new_block_size, bit_count);
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }
        // deal with the remainder in a threadblock
        if (num_remainder_in_tb > 0)
        {
            for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
            {
                prior = (data[i]) * inver_bound;
                outSize_perthread += sizeof(unsigned int);
                max = 0;
                for (j = 0; j < num_remainder_in_tb - 1; j++)
                {
                    current = (data[i + j + 1]) * inver_bound;
                    nbSampleEle++;
                    diff = current - prior;
                    prior = current;
                    if (diff > 0)
                    {
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                    else if (diff < 0)
                    {
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                }
                if (max == 0) // can be memset
                {
                    // only count the fixed length
                    outSize_perthread++;
                }
                else
                {
                    // count the fixed length, sign bits, and the saved bits
                    bit_count = (int)(log2f(max)) + 1;
                    outSize_perthread++;

                    signbytelength = estimate_signbytelength(num_remainder_in_tb - 1);
                    outSize_perthread += signbytelength;

                    savedbitsbytelength = estimate_savedbitsbytelength(num_remainder_in_tb - 1, bit_count);
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }

        // deal with the remainder outside of the normal threadblock
        if (tid == nbThreads - 1 && remainder != 0)
        {
            unsigned int num_full_block_in_rm = (remainder) / block_size;
            unsigned int num_remainder_in_rm = (remainder) % block_size;
            if (num_full_block_in_rm > 0)
            {
                // lo + 1 becuase the first one has been saved
                for (i = hi; i < nbEle - num_remainder_in_rm; i = i + block_size)
                {
                    prior = (data[i]) * inver_bound;
                    outSize_perthread += sizeof(unsigned int);
                    max = 0;
                    for (j = 0; j < new_block_size; j++)
                    {
                        current = (data[i + j + 1]) * inver_bound;
                        nbSampleEle++;
                        diff = current - prior;
                        prior = current;
                        if (diff > 0)
                        {
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                        else if (diff < 0)
                        {
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                    }
                    if (max == 0) // can be memset
                    {
                        // only count the fixed length
                        outSize_perthread++;
                    }
                    else
                    {
                        // save the fixed length, sign bits, and the saved bits
                        bit_count = (int)(log2f(max)) + 1;
                        outSize_perthread++;

                        signbytelength = estimate_signbytelength(new_block_size);
                        outSize_perthread += signbytelength;

                        savedbitsbytelength = estimate_savedbitsbytelength(new_block_size, bit_count);
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
            if (num_remainder_in_rm > 0)
            {
                // deal with the remainder in a threadblock
                for (i = nbEle - num_remainder_in_rm; i < nbEle; i = i + block_size)
                {
                    max = 0;
                    prior = (data[i]) * inver_bound;
                    outSize_perthread += sizeof(unsigned int);
                    for (j = 0; j < num_remainder_in_rm - 1; j++)
                    {
                        current = (data[i + j + 1]) * inver_bound;
                        nbSampleEle++;
                        diff = current - prior;
                        prior = current;
                        if (diff > 0)
                        {
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                        else if (diff < 0)
                        {
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                    }
                    if (max == 0) // can be memset
                    {
                        outSize_perthread++;
                    }
                    else
                    {
                        // save the fixed length, sign bits, and the saved bits
                        bit_count = (int)(log2f(max)) + 1;
                        outSize_perthread++;

                        signbytelength = estimate_signbytelength(num_remainder_in_tb - 1);
                        outSize_perthread += signbytelength;

                        savedbitsbytelength = estimate_savedbitsbytelength(num_remainder_in_tb - 1, bit_count);
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
        }

        outSize_perthread_arr[tid] = outSize_perthread;
        #pragma omp barrier
        #pragma omp single
        {
            offsets_perthread_arr[0] = 0;
            for (i = 1; i < nbThreads; i++)
            {
                offsets_perthread_arr[i] = offsets_perthread_arr[i - 1] + outSize_perthread_arr[i - 1];
                // printf("In temp arr: The offsets for thread %d is %zu bytes\n", i, offsets_perthread_arr[i]);
            }
            outSize = offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1]; 
        }
        #pragma omp barrier
        #pragma omp single
        {
            // free the buffers allocated by one thread
            free(outSize_perthread_arr);
            free(offsets_perthread_arr);
        }
    }

    // compute and return based on outSize
    size_t oriSize = sizeof(float) * nbSampleEle;
    double cr = oriSize*1.0 / outSize;
    return cr;
}

// SSIM
size_t downsampleIndex(size_t curridx, const size_t dims[], size_t n_dim, size_t sample_stride) {
    // 1D: Downsample simply by dividing the index by sample_stride.
    if (n_dim == 1) {
        return curridx / sample_stride;
    }
    // 2D: 
    // Assume dims[0] is the width and dims[1] is the height.
    // The downsample factor per axis is sqrt(sample_stride).
    else if (n_dim == 2) {
        size_t width = dims[0];
        size_t x = curridx % width;
        size_t y = curridx / width;
        // Convert sample_stride to a per–axis factor.
        size_t factor = sqrt(sample_stride);
        // Compute downsampled coordinates using integer division.
        size_t new_x = x / factor;
        size_t new_y = y / factor;
        // Calculate new width for the downsampled array.
        size_t new_width = width / factor;
        return new_y * new_width + new_x;
    }
    // 3D:
    // Assume dims[0] is width, dims[1] is height, and dims[2] is depth.
    // The per–axis factor is the cube root of sample_stride.
    else if (n_dim == 3) {
        size_t width  = dims[0];
        size_t height = dims[1];
        size_t x = curridx % width;
        size_t y = (curridx / width) % height;
        size_t z = curridx / (width * height);
        size_t factor = cbrt(sample_stride);
        size_t new_x = x / factor;
        size_t new_y = y / factor;
        size_t new_z = z / factor;
        size_t new_width  = width  / factor;
        size_t new_height = height / factor;
        return new_z * (new_width * new_height) + new_y * new_width + new_x;
    }
    return curridx;
}

struct FloatPtrArray SZp_float_estimate_ssim(float* data, size_t dims[3], size_t nbEle, size_t n_dim, int blockSize, float absErrBound, size_t sample_stride)
{
    unsigned int nbThreads = 0;
    double inver_bound = 0;
    double recipPrecision = 0;
    unsigned int threadblocksize = 0;
    unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int new_block_size = block_size - 1;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
    float eb = absErrBound;
    float* sampleData;
    float* decData;
    size_t nbSampleEle = nbEle / sample_stride;

    // size_t sample_dims[3];
    // float sample_rate = 1.0 / sample_stride;
    // for(size_t i = 0; i < 3; i++){
    //     if(n_dim == 1){
    //         sample_dims[i] = dims[i] / sample_stride;
    //     }
    //     else if(n_dim == 2){
    //         sample_dims[i] = sqrt(sample_rate)*dims[i];
    //     }
    //     else {
    //         sample_dims[i] = cbrt(sample_rate)*dims[i];
    //     }
        
    // }

    unsigned int outSize = 0;

    #pragma omp parallel
    {
        #pragma omp single
        {
            nbThreads = omp_get_num_threads();
            inver_bound = 1 / absErrBound;
            recipPrecision = 0.5f / absErrBound;
            threadblocksize = nbEle / nbThreads;
            remainder = nbEle % nbThreads;
            num_full_block_in_tb = (threadblocksize) / block_size; // minus 1 because the first one will be direcrly saved
            num_remainder_in_tb = (threadblocksize) % block_size;
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            
            // sampleData = (float *)malloc(nbSampleEle * sizeof(float));
            // decData = (float *)malloc(nbSampleEle * sizeof(float));
            // memset(sampleData, 0, nbSampleEle);
            // memset(decData, 0, nbSampleEle);
            decData = (float *)malloc(nbEle * sizeof(float));
            memset(decData, 0, nbEle);
            printf("allocsampledata\n");
        }

        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        size_t outSize_perthread = 0;
        // qutization centered with 0
        int tid = omp_get_thread_num();
        int lo = tid * threadblocksize;
        int hi = (tid + 1) * threadblocksize;

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned int signbytelength = 0; // save the bytelength
        unsigned int savedbitsbytelength = 0;

        if (num_full_block_in_tb > 0)
        {
            // iterate blocks by sample stride
            for (i = lo; i < hi - num_remainder_in_tb; i = i + block_size)//(block_size * sample_stride))
            {
                max = 0;
                prior = (data[i]) * inver_bound;
                outSize_perthread += sizeof(unsigned int);

                for (j = 0; j < new_block_size; j++)
                {
                    float curr_data = data[i + j + 1];
                    current = (data[i + j + 1]) * inver_bound;
                    
                    int currQuant = quantization_f32(curr_data, recipPrecision);
                    float pointEstimate = currQuant * 2 * eb;
                    //size_t sample_idx = i / (32 * sample_stride) * 32 + (i % 32);
                    size_t curridx = i + j + 1;
                    decData[curridx] = pointEstimate;
                    // size_t sampleidx = curridx / (block_size * sample_stride) * block_size + (curridx % block_size);
                    // // size_t sampleidx = downsampleIndex(curridx, dims, n_dim, sample_stride);
                    // if(sampleidx < nbSampleEle) {
                    //     sampleData[sampleidx] = curr_data;
                    //     decData[sampleidx] = pointEstimate;
                    // }


                    diff = current - prior;
                    prior = current;
                    
                    if (diff > 0)
                    {
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                    else if (diff < 0)
                    {
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                }
                if (max == 0) // can be memset
                {
                    // only count the fixed length
                    //sqerr 0 for constant block
                    outSize_perthread++;
                }
                else
                {
                    // count the fixed length, sign bits, and the saved bits
                    bit_count = (int)(log2f(max)) + 1;
                    outSize_perthread++;

                    signbytelength = estimate_signbytelength(new_block_size);
                    outSize_perthread += signbytelength;

                    savedbitsbytelength = estimate_savedbitsbytelength(new_block_size, bit_count);
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }
        printf("startremainder\n");
        // deal with the remainder in a threadblock
        if (false && num_remainder_in_tb > 0)
        {
            for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
            {
                prior = (data[i]) * inver_bound;
                outSize_perthread += sizeof(unsigned int);
                max = 0;
                for (j = 0; j < num_remainder_in_tb - 1; j++)
                {
                    float curr_data = data[i + j + 1];
                    current = (data[i + j + 1]) * inver_bound;

                    int currQuant = quantization_f32(curr_data, recipPrecision);
                    float pointEstimate = currQuant * 2 * eb;

                    size_t curridx = i + j + 1;
                    decData[curridx] = pointEstimate;
                    // size_t sampleidx = curridx / (block_size * sample_stride) * block_size + (curridx % block_size);
                    // // size_t sampleidx = downsampleIndex(curridx, dims, n_dim, sample_stride);
                    // if(sampleidx < nbSampleEle) {
                    //     sampleData[sampleidx] = curr_data;
                    //     decData[sampleidx] = pointEstimate;
                    // }

                    // size_t curridx = i + j + 1;
                    // size_t sampleidx = curridx + (curridx % block_size);
                    // if(sampleidx < nbSampleEle) {
                    //     sampleData[sampleidx] = curr_data;
                    //     decData[sampleidx] = pointEstimate;
                    // }

                    diff = current - prior;
                    prior = current;

                    if (diff > 0)
                    {
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                    else if (diff < 0)
                    {
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                }
                if (max == 0) // can be memset
                {
                    // only count the fixed length
                    outSize_perthread++;
                }
                else
                {
                    // count the fixed length, sign bits, and the saved bits
                    bit_count = (int)(log2f(max)) + 1;
                    outSize_perthread++;

                    signbytelength = estimate_signbytelength(num_remainder_in_tb - 1);
                    outSize_perthread += signbytelength;

                    savedbitsbytelength = estimate_savedbitsbytelength(num_remainder_in_tb - 1, bit_count);
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }
        printf("outthread\n");
        // deal with the remainder outside of the normal threadblock
        if (false && tid == nbThreads - 1 && remainder != 0)
        {
            unsigned int num_full_block_in_rm = (remainder) / block_size;
            unsigned int num_remainder_in_rm = (remainder) % block_size;
            if (num_full_block_in_rm > 0)
            {
                // lo + 1 becuase the first one has been saved
                for (i = hi; i < nbEle - num_remainder_in_rm; i = i + block_size)
                {
                    prior = (data[i]) * inver_bound;
                    outSize_perthread += sizeof(unsigned int);
                    max = 0;
                    for (j = 0; j < new_block_size; j++)
                    {
                        float curr_data = data[i + j + 1];
                        current = (data[i + j + 1]) * inver_bound;

                        int currQuant = quantization_f32(curr_data, recipPrecision);
                        float pointEstimate = currQuant * 2 * eb;

                        size_t curridx = i + j + 1;
                        decData[curridx] = pointEstimate;
                        // size_t sampleidx = curridx / (block_size * sample_stride) * block_size + (curridx % block_size);
                        // // size_t sampleidx = downsampleIndex(curridx, dims, n_dim, sample_stride);
                        // if(sampleidx < nbSampleEle) {
                        //     sampleData[sampleidx] = curr_data;
                        //     decData[sampleidx] = pointEstimate;
                        // }
                        
                        diff = current - prior;
                        prior = current;

                        if (diff > 0)
                        {
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                        else if (diff < 0)
                        {
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                    }
                    if (max == 0) // can be memset
                    {
                        // only count the fixed length
                        outSize_perthread++;
                    }
                    else
                    {
                        // save the fixed length, sign bits, and the saved bits
                        bit_count = (int)(log2f(max)) + 1;
                        outSize_perthread++;

                        signbytelength = estimate_signbytelength(new_block_size);
                        outSize_perthread += signbytelength;

                        savedbitsbytelength = estimate_savedbitsbytelength(new_block_size, bit_count);
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
            if (num_remainder_in_rm > 0)
            {
                // deal with the remainder in a threadblock
                for (i = nbEle - num_remainder_in_rm; i < nbEle; i = i + block_size)
                {
                    max = 0;
                    prior = (data[i]) * inver_bound;
                    outSize_perthread += sizeof(unsigned int);
                    for (j = 0; j < num_remainder_in_rm - 1; j++)
                    {
                        float curr_data = data[i + j + 1];
                        current = (data[i + j + 1]) * inver_bound;

                        int currQuant = quantization_f32(curr_data, recipPrecision);
                        float pointEstimate = currQuant * 2 * eb;

                        size_t curridx = i + j + 1;
                        decData[curridx] = pointEstimate;
                        // size_t sampleidx = curridx / (block_size * sample_stride) * block_size + (curridx % block_size);
                        // // size_t sampleidx = downsampleIndex(curridx, dims, n_dim, sample_stride);
                        // if(sampleidx < nbSampleEle) {
                        //     sampleData[sampleidx] = curr_data;
                        //     decData[sampleidx] = pointEstimate;
                        // }
                        
                        diff = current - prior;
                        prior = current;

                        if (diff > 0)
                        {
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                        else if (diff < 0)
                        {
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                    }
                    if (max == 0) // can be memset
                    {
                        outSize_perthread++;
                    }
                    else
                    {
                        // save the fixed length, sign bits, and the saved bits
                        bit_count = (int)(log2f(max)) + 1;
                        outSize_perthread++;

                        signbytelength = estimate_signbytelength(num_remainder_in_tb - 1);
                        outSize_perthread += signbytelength;

                        savedbitsbytelength = estimate_savedbitsbytelength(num_remainder_in_tb - 1, bit_count);
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
        }

        outSize_perthread_arr[tid] = outSize_perthread;

        // #pragma omp barrier
        // #pragma omp single
        // {
        //     offsets_perthread_arr[0] = 0;
        //     for (i = 1; i < nbThreads; i++)
        //     {
        //         offsets_perthread_arr[i] = offsets_perthread_arr[i - 1] + outSize_perthread_arr[i - 1];
        //         // printf("In temp arr: The offsets for thread %d is %zu bytes\n", i, offsets_perthread_arr[i]);
        //     }
        //     outSize = offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1]; 

        // }
        #pragma omp barrier
        #pragma omp single
        {
            // free the buffers allocated by one thread
            free(outSize_perthread_arr);
            free(offsets_perthread_arr);
        }
    }
    printf("resultgive\n");

    struct FloatPtrArray result;
    result.arr[0] = sampleData;
    result.arr[1] = decData;
    
    return result;
}

// PSNR
double SZp_float_estimate_psnr(float* data, size_t nbEle, int blockSize, float absErrBound, size_t sample_stride)
{
    unsigned int nbThreads = 0;
    double inver_bound = 0;
    double recipPrecision = 0;
    unsigned int threadblocksize = 0;
    unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int new_block_size = block_size - 1;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
    double* sqerr_perthread_arr;
    float* max_perthread_arr;
    float* min_perthread_arr;
    float eb = absErrBound;

    size_t nbSampleEle = 0;

    unsigned int outSize = 0;
    double global_sqerr;
    float global_max;
    float global_min;

    #pragma omp parallel
    {
        #pragma omp single
        {
            nbThreads = omp_get_num_threads();
            inver_bound = 1 / absErrBound;
            recipPrecision = 0.5f / absErrBound;
            threadblocksize = nbEle / nbThreads;
            remainder = nbEle % nbThreads;
            num_full_block_in_tb = (threadblocksize) / block_size; // minus 1 because the first one will be direcrly saved
            num_remainder_in_tb = (threadblocksize) % block_size;
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            
            sqerr_perthread_arr = (double *)malloc(nbThreads * sizeof(double));
            max_perthread_arr = (float *)malloc(nbThreads * sizeof(float));
            min_perthread_arr = (float *)malloc(nbThreads * sizeof(float));
        }

        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        size_t outSize_perthread = 0;
        // qutization centered with 0
        int tid = omp_get_thread_num();
        int lo = tid * threadblocksize;
        int hi = (tid + 1) * threadblocksize;

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned int signbytelength = 0; // save the bytelength
        unsigned int savedbitsbytelength = 0;

        // init perthread vars
        double local_sqerr = 0;
        float local_max = FLT_MIN;
        float local_min = FLT_MAX;

        if (num_full_block_in_tb > 0)
        {
            // iterate blocks by sample stride
            for (i = lo; i < hi - num_remainder_in_tb; i = i + (block_size * sample_stride))
            {
                max = 0;
                prior = (data[i]) * inver_bound;
                outSize_perthread += sizeof(unsigned int);

                for (j = 0; j < new_block_size; j++)
                {
                    float curr_data = data[i + j + 1];
                    current = (data[i + j + 1]) * inver_bound;
                    
                    int currQuant = quantization_f32(curr_data, recipPrecision);
                    float pointEstimate = currQuant * 2 * eb;
                    local_sqerr += pow(curr_data - pointEstimate, 2); 
                    nbSampleEle++;

                    diff = current - prior;
                    prior = current;
                    if (curr_data < local_min)
                    {
                        local_min = curr_data;
                    }
                    if (curr_data > local_max)
                    {
                        local_max = curr_data;
                    }
                    
                    if (diff > 0)
                    {
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                    else if (diff < 0)
                    {
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                }
                if (max == 0) // can be memset
                {
                    // only count the fixed length
                    //sqerr 0 for constant block
                    outSize_perthread++;
                }
                else
                {
                    // count the fixed length, sign bits, and the saved bits
                    bit_count = (int)(log2f(max)) + 1;
                    outSize_perthread++;

                    signbytelength = estimate_signbytelength(new_block_size);
                    outSize_perthread += signbytelength;

                    savedbitsbytelength = estimate_savedbitsbytelength(new_block_size, bit_count);
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }
        // deal with the remainder in a threadblock
        if (num_remainder_in_tb > 0)
        {
            for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
            {
                prior = (data[i]) * inver_bound;
                outSize_perthread += sizeof(unsigned int);
                max = 0;
                for (j = 0; j < num_remainder_in_tb - 1; j++)
                {
                    float curr_data = data[i + j + 1];
                    current = (data[i + j + 1]) * inver_bound;

                    int currQuant = quantization_f32(curr_data, recipPrecision);
                    float pointEstimate = currQuant * 2 * eb;
                    local_sqerr += pow(curr_data - pointEstimate, 2);

                    nbSampleEle++;
                    diff = current - prior;
                    prior = current;
                    if (curr_data < local_min)
                    {
                        local_min = curr_data;
                    }
                    if (curr_data > local_max)
                    {
                        local_max = curr_data;
                    }

                    if (diff > 0)
                    {
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                    else if (diff < 0)
                    {
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                    }
                }
                if (max == 0) // can be memset
                {
                    // only count the fixed length
                    outSize_perthread++;
                }
                else
                {
                    // count the fixed length, sign bits, and the saved bits
                    bit_count = (int)(log2f(max)) + 1;
                    outSize_perthread++;

                    signbytelength = estimate_signbytelength(num_remainder_in_tb - 1);
                    outSize_perthread += signbytelength;

                    savedbitsbytelength = estimate_savedbitsbytelength(num_remainder_in_tb - 1, bit_count);
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }

        // deal with the remainder outside of the normal threadblock
        if (tid == nbThreads - 1 && remainder != 0)
        {
            unsigned int num_full_block_in_rm = (remainder) / block_size;
            unsigned int num_remainder_in_rm = (remainder) % block_size;
            if (num_full_block_in_rm > 0)
            {
                // lo + 1 becuase the first one has been saved
                for (i = hi; i < nbEle - num_remainder_in_rm; i = i + block_size)
                {
                    prior = (data[i]) * inver_bound;
                    outSize_perthread += sizeof(unsigned int);
                    max = 0;
                    for (j = 0; j < new_block_size; j++)
                    {
                        float curr_data = data[i + j + 1];
                        current = (data[i + j + 1]) * inver_bound;

                        int currQuant = quantization_f32(curr_data, recipPrecision);
                        float pointEstimate = currQuant * 2 * eb;
                        local_sqerr += pow(curr_data - pointEstimate, 2); 

                        nbSampleEle++;
                        diff = current - prior;
                        prior = current;
                        if (curr_data < local_min)
                        {
                            local_min = curr_data;
                        }
                        if (curr_data > local_max)
                        {
                            local_max = curr_data;
                        }

                        if (diff > 0)
                        {
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                        else if (diff < 0)
                        {
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                    }
                    if (max == 0) // can be memset
                    {
                        // only count the fixed length
                        outSize_perthread++;
                    }
                    else
                    {
                        // save the fixed length, sign bits, and the saved bits
                        bit_count = (int)(log2f(max)) + 1;
                        outSize_perthread++;

                        signbytelength = estimate_signbytelength(new_block_size);
                        outSize_perthread += signbytelength;

                        savedbitsbytelength = estimate_savedbitsbytelength(new_block_size, bit_count);
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
            if (num_remainder_in_rm > 0)
            {
                // deal with the remainder in a threadblock
                for (i = nbEle - num_remainder_in_rm; i < nbEle; i = i + block_size)
                {
                    max = 0;
                    prior = (data[i]) * inver_bound;
                    outSize_perthread += sizeof(unsigned int);
                    for (j = 0; j < num_remainder_in_rm - 1; j++)
                    {
                        float curr_data = data[i + j + 1];
                        current = (data[i + j + 1]) * inver_bound;

                        int currQuant = quantization_f32(curr_data, recipPrecision);
                        float pointEstimate = currQuant * 2 * eb;
                        local_sqerr += pow(curr_data - pointEstimate, 2); 

                        nbSampleEle++;
                        diff = current - prior;
                        prior = current;
                        if (curr_data < local_min)
                        {
                            local_min = curr_data;
                        }
                        if (curr_data > local_max)
                        {
                            local_max = curr_data;
                        }

                        if (diff > 0)
                        {
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                        else if (diff < 0)
                        {
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                        }
                    }
                    if (max == 0) // can be memset
                    {
                        outSize_perthread++;
                    }
                    else
                    {
                        // save the fixed length, sign bits, and the saved bits
                        bit_count = (int)(log2f(max)) + 1;
                        outSize_perthread++;

                        signbytelength = estimate_signbytelength(num_remainder_in_tb - 1);
                        outSize_perthread += signbytelength;

                        savedbitsbytelength = estimate_savedbitsbytelength(num_remainder_in_tb - 1, bit_count);
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
        }

        outSize_perthread_arr[tid] = outSize_perthread;
        sqerr_perthread_arr[tid] = local_sqerr;
        max_perthread_arr[tid] = local_max;
        min_perthread_arr[tid] = local_min;

        #pragma omp barrier
        #pragma omp single
        {
            offsets_perthread_arr[0] = 0;
            for (i = 1; i < nbThreads; i++)
            {
                offsets_perthread_arr[i] = offsets_perthread_arr[i - 1] + outSize_perthread_arr[i - 1];
                // printf("In temp arr: The offsets for thread %d is %zu bytes\n", i, offsets_perthread_arr[i]);
            }
            outSize = offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1]; 

            global_sqerr = sqerr_perthread_arr[0];
            global_max = max_perthread_arr[0];
            global_min = min_perthread_arr[0];
            for(int z = 1; z < nbThreads; z++) {
                global_sqerr += sqerr_perthread_arr[z];
                if(max_perthread_arr[z] > global_max)
                    global_max = max_perthread_arr[z];
                if(min_perthread_arr[z] > global_min)
                    global_min = min_perthread_arr[z];
            }
        }
        #pragma omp barrier
        #pragma omp single
        {
            // free the buffers allocated by one thread
            free(outSize_perthread_arr);
            free(offsets_perthread_arr);

            free(sqerr_perthread_arr);
            free(max_perthread_arr);
            free(min_perthread_arr);
        }
    }

    // compute psnr
    double eps = 1e-16;
    double mse = global_sqerr / nbSampleEle;
    double value_range = global_max - global_min;
    double psnr = -20.0*log10((sqrt(mse)/value_range) + eps);
    
    return psnr;
}
