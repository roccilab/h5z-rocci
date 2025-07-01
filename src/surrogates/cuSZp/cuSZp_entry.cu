#include <cuSZp/cuSZp_f32.h>
#include <iostream>
#include "qcat_ssim.h"

__global__ void estimate_SZp_compress_kernel_f32(const float* const __restrict__ oriData, double* total_compressed_size, const float eb, const size_t nbEle, const size_t sample_stride);
__global__ void estimate_SZp_psnr_kernel_f32(const float* const __restrict__ oriData, float* maxVals, float* minVals, double* square_errors, const float eb, const size_t nbEle, const size_t sample_stride);
__global__ void estimate_SZp_ssim_kernel_f32(const float* const __restrict__ oriData, size_t n_dim, size_t dims[3], float* sampleData, float* decData, const float eb, const size_t nbEle, const size_t sample_stride);


double SZp_estimate_compress_hostptr_f32(float* oriData, size_t nbEle, float errorBound, size_t sample_stride)
{
    // Data blocking.
    int bsize = cmp_tblock_size_f32;
    int gsize = (nbEle + bsize * cmp_chunk_f32 - 1) / (bsize * cmp_chunk_f32);
    int sample_gsize = ((nbEle / sample_stride) + bsize * cmp_chunk_f32 - 1) / (bsize * cmp_chunk_f32);
    int cmpOffSize = gsize + 1;
    int pad_nbEle = gsize * bsize * cmp_chunk_f32;
    if (pad_nbEle > nbEle){
        pad_nbEle = nbEle;
    }

    size_t BLOCKSIZE = 32;
    size_t nbBlocks = pad_nbEle / BLOCKSIZE;
    size_t nbChunks = pad_nbEle / cmp_chunk_f32;

    // Initializing global memory for GPU compression.
    float* d_oriData;
    double* d_total_compressed_size;
    cudaMalloc((void**)&d_oriData, sizeof(float)*pad_nbEle);
    cudaMemcpy(d_oriData, oriData, sizeof(float)*pad_nbEle, cudaMemcpyHostToDevice);
    // cudaMalloc(&d_total_compressed_size, sizeof(double) * nbBlocks);
    // cudaMemset(d_total_compressed_size, 0, sizeof(double) * nbBlocks);
    cudaMalloc(&d_total_compressed_size, sizeof(double) * nbChunks);
    cudaMemset(d_total_compressed_size, 0, sizeof(double) * nbChunks);
    cudaMemset(d_oriData + nbEle, 0, (pad_nbEle - nbEle) * sizeof(float));

    // Initializing CUDA Stream.
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // cuSZp GPU compression.
    dim3 blockSize(bsize);
    dim3 gridSize(sample_gsize); // only launch sample grid num blocks
    estimate_SZp_compress_kernel_f32<<<gridSize, blockSize, 0, stream>>>(d_oriData, d_total_compressed_size, errorBound, nbEle, sample_stride);
    cudaDeviceSynchronize();

    // Obtain compression ratio and move data back to CPU. 
    // double* h_total_compressed_size = (double*) malloc(sizeof(double) * nbBlocks);
    // cudaMemcpy(h_total_compressed_size, d_total_compressed_size, sizeof(double)*nbBlocks, cudaMemcpyDeviceToHost); 
    double* h_total_compressed_size = (double*) malloc(sizeof(double) * nbChunks);
    cudaMemcpy(h_total_compressed_size, d_total_compressed_size, sizeof(double)*nbChunks, cudaMemcpyDeviceToHost); 


    // compute on chunk stride instead of block stride
    // size_t nbSampleBlocks = (nbBlocks / sample_stride);
    // size_t nbSampleEle = nbSampleBlocks * BLOCKSIZE;
    // double original_size = nbSampleEle * sizeof(float);
    // size_t nbChunk = nbEle / cmp_chunk_f32;
    size_t nbSampleChunk = nbChunks / sample_stride;
    size_t nbSampleEle = nbSampleChunk * cmp_chunk_f32;
    double original_size = nbSampleEle * sizeof(float);

    // sum chunk compressed sizes
    double total_size = 0;
    // for(int i = 0; i < nbBlocks; i+=sample_stride) {
    for(int i = 0; i < nbChunks; i+=sample_stride) {
        total_size += h_total_compressed_size[i];
    }

    double compression_ratio = original_size / (total_size);

    // Free memory that is used.
    cudaFree(d_oriData);
    cudaFree(d_total_compressed_size);
    free(h_total_compressed_size);
    cudaStreamDestroy(stream);

    return compression_ratio;
}


double SZp_estimate_psnr_hostptr_f32(float* oriData, size_t nbEle, float errorBound, size_t sample_stride)
{
    // Data blocking.
    int bsize = cmp_tblock_size_f32;
    int gsize = (nbEle + bsize * cmp_chunk_f32 - 1) / (bsize * cmp_chunk_f32);
    int cmpOffSize = gsize + 1;
    int pad_nbEle = gsize * bsize * cmp_chunk_f32;
    if (pad_nbEle > nbEle){
        pad_nbEle = nbEle;
    }
    int num_chunks = pad_nbEle / cmp_chunk_f32;

    size_t BLOCKSIZE = 32;
    size_t nbBlocks = pad_nbEle / BLOCKSIZE;

    // Initializing global memory for GPU compression.
    float* d_oriData;
    float* d_maxVals;
    float* d_minVals;
    double* d_errors;
    cudaMalloc((void**)&d_oriData, sizeof(float)*pad_nbEle);
    cudaMemcpy(d_oriData, oriData, sizeof(float)*pad_nbEle, cudaMemcpyHostToDevice);
    cudaMalloc(&d_maxVals, sizeof(float) * num_chunks);
    cudaMemset(d_maxVals, 0, sizeof(float) * num_chunks);
    cudaMalloc(&d_minVals, sizeof(float) * num_chunks);
    cudaMemset(d_minVals, 0, sizeof(float) * num_chunks);
    cudaMalloc(&d_errors, sizeof(double) * num_chunks);
    cudaMemset(d_errors, 0, sizeof(double) * num_chunks);
    cudaMemset(d_oriData + nbEle, 0, (pad_nbEle - nbEle) * sizeof(float));

    // Initializing CUDA Stream.
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // cuSZp GPU compression.
    dim3 blockSize(bsize);
    dim3 gridSize(gsize);
    estimate_SZp_psnr_kernel_f32<<<gridSize, blockSize, 0, stream>>>(d_oriData, d_maxVals, d_minVals, d_errors, errorBound, nbEle, sample_stride);
    cudaDeviceSynchronize();

    // Obtain compression ratio and move data back to CPU. 
    float* h_maxVals = (float*) malloc(sizeof(float) * num_chunks);
    cudaMemcpy(h_maxVals, d_maxVals, sizeof(float)*num_chunks, cudaMemcpyDeviceToHost); 
    float* h_minVals = (float*) malloc(sizeof(float) * num_chunks);
    cudaMemcpy(h_minVals, d_minVals, sizeof(float)*num_chunks, cudaMemcpyDeviceToHost); 
    double* h_errors = (double*) malloc(sizeof(double) * num_chunks);
    cudaMemcpy(h_errors, d_errors, sizeof(double)*num_chunks, cudaMemcpyDeviceToHost); 

    size_t nbSampleBlocks = (nbBlocks / sample_stride);
    size_t nbSampleEle = nbSampleBlocks * BLOCKSIZE;
    double original_size = nbSampleEle * sizeof(float);

    // aggregate values
    float maxVal = h_maxVals[0];
    float minVal = h_minVals[0];
    double sq_err = 0;
    for(int i = 0; i < num_chunks; i++) {
        if(h_maxVals[i] > maxVal) maxVal = h_maxVals[i];
        if(h_minVals[i] < minVal) minVal = h_minVals[i];
        sq_err += h_errors[i];
    }

    double valRange = maxVal - minVal;
    double mse = sq_err / nbSampleEle;
    double eps = 1e-16;
    double psnr = -20.0*log10((sqrt(mse) / valRange) + eps);

    // Free memory that is used.
    cudaFree(d_oriData);
    cudaFree(d_maxVals);
    cudaFree(d_minVals);
    cudaFree(d_errors);
    free(h_maxVals);
    free(h_minVals);
    free(h_errors);
    cudaStreamDestroy(stream);

    return psnr;
}

double SZp_estimate_ssim_hostptr_f32(float* oriData, size_t dims[3], size_t n_dim, size_t nbEle, float errorBound, size_t sample_stride)
{
    // Data blocking.
    int bsize = cmp_tblock_size_f32;
    int gsize = (nbEle + bsize * cmp_chunk_f32 - 1) / (bsize * cmp_chunk_f32);
    int cmpOffSize = gsize + 1;
    int pad_nbEle = gsize * bsize * cmp_chunk_f32;
    if (pad_nbEle > nbEle){
        pad_nbEle = nbEle;
    }

    size_t BLOCKSIZE = 32;
    size_t nbBlocks = pad_nbEle / BLOCKSIZE;
    size_t nbSampleBlocks = (nbBlocks / sample_stride);
    size_t nbSampleEle = nbSampleBlocks * BLOCKSIZE;
    float sample_rate = 1.0 / sample_stride;

    size_t sample_dims[3];
    size_t num_sample_points = 1;
    // for(size_t i = 0; i < 3; i++){
    //     if(n_dim == 1){
    //         sample_dims[i] = dims[i] / sample_stride;
    //     }
    //     else if(n_dim == 2){
    //         sample_dims[i] = sqrt(sample_rate)*dims[i];
    //     }
    //     else {
    //         sample_dims[i] = std::cbrt(sample_rate)*dims[i];
    //     }
    //     if(sample_dims[i] > 0){
    //         num_sample_points *= sample_dims[i];
    //     }
    // }
    // nbSampleEle = num_sample_points;

    if(n_dim == 1){
        sample_dims[0] = dims[0] / (float) sample_stride;
    }
    else if(n_dim == 2){
        sample_dims[0] = dims[0] / (float) sample_stride;
        sample_dims[1] = dims[1];
    }
    else {
        sample_dims[0] = dims[0] / (float) sample_stride;
        sample_dims[1] = dims[1];
        sample_dims[2] = dims[2];
    }

    printf("Ndim %i, SS: %i, dims: %i %i %i\n", n_dim, sample_stride, sample_dims[0], sample_dims[1], sample_dims[2]);

    // Initializing global memory for GPU compression.
    float* d_oriData;
    float* d_sampleData;
    float* d_decData;
    cudaMalloc((void**)&d_oriData, sizeof(float)*pad_nbEle);
    cudaMemcpy(d_oriData, oriData, sizeof(float)*pad_nbEle, cudaMemcpyHostToDevice);
    cudaMalloc(&d_sampleData, sizeof(float) * nbSampleEle);
    cudaMemset(d_sampleData, 0, sizeof(float) * nbSampleEle);
    cudaMalloc(&d_decData, sizeof(float) * nbSampleEle);
    cudaMemset(d_decData, 0, sizeof(float) * nbSampleEle);
    cudaMemset(d_oriData + nbEle, 0, (pad_nbEle - nbEle) * sizeof(float));

    // Initializing CUDA Stream.
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // cuSZp GPU compression.
    dim3 blockSize(bsize);
    dim3 gridSize(gsize);
    estimate_SZp_ssim_kernel_f32<<<gridSize, blockSize, 0, stream>>>(d_oriData, n_dim, dims, d_sampleData, d_decData, errorBound, nbEle, sample_stride);
    cudaDeviceSynchronize();

    // Obtain compression ratio and move data back to CPU. 
    float* h_sampleData = (float*) malloc(sizeof(float) * nbSampleEle);
    cudaMemcpy(h_sampleData, d_sampleData, sizeof(float)*nbSampleEle, cudaMemcpyDeviceToHost); 
    float* h_decData = (float*) malloc(sizeof(float) * nbSampleEle);
    cudaMemcpy(h_decData, d_decData, sizeof(float)*nbSampleEle, cudaMemcpyDeviceToHost); 

    double ssim;
	if(n_dim == 1){
		ssim = SSIM_1d_windowed_float(h_sampleData, h_decData, sample_dims[0], 8, 8);
	}
	else if (n_dim == 2) {
		ssim = SSIM_2d_windowed_float(h_sampleData, h_decData, sample_dims[1], sample_dims[0], 8,8, 8,8);
	}
	else {
		ssim = SSIM_3d_windowed_float(h_sampleData, h_decData, sample_dims[2], sample_dims[1], sample_dims[0], 8,8,8, 8,8,8);
	}

    // Free memory that is used.
    cudaFree(d_oriData);
    cudaFree(d_sampleData);
    cudaFree(d_decData);
    free(h_sampleData);
    free(h_decData);
    cudaStreamDestroy(stream);

    return ssim;
}
