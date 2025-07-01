#include <cuSZp/cuSZp_f32.h>

__device__ inline int quantization_f32(float data, float recipPrecision)
{
    float dataRecip = data*recipPrecision;
    int s = dataRecip>=-0.5f?0:1;
    return (int)(dataRecip+0.5f) - s;
}


__device__ inline int get_bit_num(unsigned int x)
{
    return (sizeof(unsigned int)*8) - __clz(x);
}

// CR
__global__ void estimate_SZp_compress_kernel_f32(const float* const __restrict__ oriData, double* total_compressed_size, const float eb, const size_t nbEle, const size_t sample_stride)
{
    __shared__ unsigned int base_idx;

    const int tid = threadIdx.x;
    // scale by stride
    const int idx = (blockIdx.x * blockDim.x + tid) * sample_stride;
    const int lane = idx & 31;
    const int warp = idx >> 5;
    const int block_num = cmp_chunk_f32/32;
    const int start_idx = idx * cmp_chunk_f32;
    const int start_block_idx = start_idx/32;
    const int rate_ofs = (nbEle+31)/32;
    const float recipPrecision = 0.5f/eb;

    int temp_start_idx, temp_end_idx;
    int quant_chunk_idx;
    int block_idx;
    int currQuant, lorenQuant, prevQuant, maxQuant;
    int absQuant[cmp_chunk_f32];
    unsigned int thread_ofs = 0;

    size_t num_zero_blocks = 0;
    double local_compressed_size = 0;

    // changes to chunk stride
    if(start_idx + cmp_chunk_f32 < nbEle){

        // iterate by sample stride to skip blocks
        // for(int j=0; j<block_num; j+=sample_stride)
        for(int j=0; j<block_num; j++)
        {
            temp_start_idx = start_idx + j*32;
            temp_end_idx = temp_start_idx + 32;
            block_idx = start_block_idx+j;

            double compressed_size = 0;

            prevQuant = 0;
            maxQuant = 0;
            
            bool is_zero_block = true;
            bool is_zero_quant_block = true;
            for(int i=temp_start_idx; i<temp_end_idx; i++){
                // compute lorenzo quantization index
                float curr = oriData[i];
                if(curr != 0) is_zero_block = false;
                currQuant = quantization_f32(curr, recipPrecision);
                if(currQuant != 0) is_zero_quant_block = false;
                lorenQuant = currQuant - prevQuant;
                prevQuant = currQuant;
                int absQuant = abs(lorenQuant);
                // update max required bit size for this block
                maxQuant = maxQuant > absQuant ? maxQuant : absQuant;   
            }

            // compute compressed size of this block in bytes
            if(is_zero_block || is_zero_quant_block){
                num_zero_blocks++;
                compressed_size += 1;
            }
            else {
                int block_bit_size = get_bit_num(maxQuant);
                compressed_size += (block_bit_size + 1) * 32 / 8;
            }

            // write out compressed size for this block
            // total_compressed_size[block_idx] = compressed_size;
            local_compressed_size += compressed_size;

        }

        total_compressed_size[idx] = local_compressed_size;
    }

}

// PSNR
__global__ void estimate_SZp_psnr_kernel_f32(const float* const __restrict__ oriData, float* max_vals, float* min_vals, double* square_errors, const float eb, const size_t nbEle, const size_t sample_stride)
{
    __shared__ unsigned int base_idx;

    const int tid = threadIdx.x;
    const int idx = blockIdx.x * blockDim.x + tid;
    const int lane = idx & 31;
    const int warp = idx >> 5;
    const int block_num = cmp_chunk_f32/32;
    const int start_idx = idx * cmp_chunk_f32;
    const int start_block_idx = start_idx/32;
    const int rate_ofs = (nbEle+31)/32;
    const float recipPrecision = 0.5f/eb;

    int temp_start_idx, temp_end_idx;
    int quant_chunk_idx;
    int block_idx;
    int currQuant, lorenQuant, prevQuant, maxQuant;
    int absQuant[cmp_chunk_f32];
    unsigned int thread_ofs = 0;

    size_t num_zero_blocks = 0;

    double sq_err = 0;
    float minVal = oriData[start_idx];
    float maxVal = oriData[start_idx];

    // iterate by sample stride to skip blocks
    for(int j=0; j<block_num; j+=sample_stride)
    {
        temp_start_idx = start_idx + j*32;
        temp_end_idx = temp_start_idx + 32;
        block_idx = start_block_idx+j;

        double compressed_size = 0;

        prevQuant = 0;
        maxQuant = 0;
        
        bool is_zero_block = true;
        bool is_zero_quant_block = true;
        for(int i=temp_start_idx; i<temp_end_idx; i++){
            // compute lorenzo quantization index
            float curr = oriData[i];
			if(curr < minVal) minVal = curr;
			if(curr > maxVal) maxVal = curr;

            currQuant = quantization_f32(curr, recipPrecision);
			float pointEstimate = currQuant * 2 * eb;
			sq_err += pow(curr - pointEstimate, 2);   
        }

    }

    max_vals[idx] = maxVal;
    min_vals[idx] = minVal;
    square_errors[idx] = sq_err;

}

__device__ size_t get_sampleidx(size_t idx, size_t sample_stride, size_t n_dim) 
{
    size_t sample_idx = idx;
    float sample_rate = 1.0 / sample_stride;
    if(n_dim == 1){
        sample_idx = sample_idx * sample_rate;
    }
    else if(n_dim == 2){
        sample_idx = sample_idx * sqrt(sample_rate);
    }
    else {
        sample_idx = sample_idx * std::cbrt(sample_rate);
    }
    return sample_idx;
}

// SSIM
// __global__ void estimate_SZp_ssim_kernel_f32(const float* const __restrict__ oriData, size_t n_dim, float* sampleData, float* decData, const float eb, const size_t nbEle, const size_t sample_stride)
// {
//     __shared__ unsigned int base_idx;

//     const int tid = threadIdx.x;
//     const int idx = blockIdx.x * blockDim.x + tid;
//     const int lane = idx & 31;
//     const int warp = idx >> 5;
//     const int block_num = cmp_chunk_f32/32;
//     const int start_idx = idx * cmp_chunk_f32;
//     const int start_block_idx = start_idx/32;
//     const int rate_ofs = (nbEle+31)/32;
//     const float recipPrecision = 0.5f/eb;

//     size_t nbBlocks = nbEle / 32;
//     size_t nbSampleBlocks = (nbBlocks / sample_stride);
//     size_t nbSampleEle = nbSampleBlocks * 32;

//     int temp_start_idx, temp_end_idx;
//     int quant_chunk_idx;
//     int block_idx;
//     int currQuant, lorenQuant, prevQuant, maxQuant;
//     int absQuant[cmp_chunk_f32];
//     unsigned int thread_ofs = 0;

//     size_t num_zero_blocks = 0;
//     // set proper start position in sample array
//     // size_t sample_idx = start_idx;
//     // float sample_rate = 1.0 / sample_stride;
//     // if(n_dim == 1){
//     //     sample_idx = sample_idx * sample_rate;
//     // }
//     // else if(n_dim == 2){
//     //     sample_idx = sample_idx * sqrt(sample_rate);
//     // }
//     // else {
//     //     sample_idx = sample_idx * std::cbrt(sample_rate);
//     // }

//     // if(sample_idx < nbSampleEle){

//         // iterate by sample stride to skip blocks
//         for(int j=0; j<block_num; j+=sample_stride)
//         {
//             temp_start_idx = start_idx + j*32;
//             temp_end_idx = temp_start_idx + 32;
//             block_idx = start_block_idx+j;

//             size_t sample_idx = get_sampleidx(temp_start_idx, sample_stride, n_dim);

//             double compressed_size = 0;

//             prevQuant = 0;
//             maxQuant = 0;
            
//             bool is_zero_block = true;
//             bool is_zero_quant_block = true;
//             for(int i=temp_start_idx; i<temp_end_idx; i++){
                
//                 // compute lorenzo quantization index
//                 float curr = oriData[i];
//                 currQuant = quantization_f32(curr, recipPrecision);
//                 float pointEstimate = currQuant * 2 * eb;
                
//                 sampleData[sample_idx] = curr;
//                 decData[sample_idx] = pointEstimate;
//                 sample_idx++;

//                 if(sample_idx >= nbSampleEle) break;
//             }

//             // if(sample_idx >= nbSampleEle) break;

//         }
//     // }

// }

// __global__ void estimate_SZp_ssim_kernel_f32(const float* oriData, size_t n_dim, float* sampleData, float* decData, 
//                                              const float eb, const size_t nbEle, const size_t sample_stride)
// {
//     const int tid = blockIdx.x * blockDim.x + threadIdx.x;
//     const float recipPrecision = 0.5f / eb;
//     size_t nbBlocks = nbEle / 32;
//     size_t nbSampleBlocks = (nbBlocks / sample_stride);
//     size_t nbSampleEle = nbSampleBlocks * 32;

//     for (size_t i = tid; i < nbEle; i += gridDim.x * blockDim.x) {
//         if (i % (32 * sample_stride) < 32) {  // Only process every sample_stride'th block
//             size_t sample_idx = i / (32 * sample_stride) * 32 + (i % 32);
            
//             if (sample_idx < nbSampleEle) {
//                 float curr = oriData[i];
//                 int currQuant = quantization_f32(curr, recipPrecision);
//                 float pointEstimate = currQuant * 2 * eb;
                
//                 sampleData[sample_idx] = curr;
//                 decData[sample_idx] = pointEstimate;
//             }
//         }
//     }
// }

__device__ size_t multidim_sample_idx(size_t idx, size_t n_dim, size_t dims[3], size_t sample_dims[3], size_t sample_stride) {
    float sample_rate = 1 / sample_stride;
    if (n_dim == 1){
        return idx / sample_stride;
    }
    else if(n_dim == 2){
        size_t x = idx % dims[0];
        size_t y = idx / dims[1];
        size_t sample_x = x / (1/sqrt(sample_rate));
        size_t sample_y = y / (1/sqrt(sample_rate));
        return sample_x + sample_dims[0] * sample_y;
    }
    else {
        size_t x = idx % dims[0];
        size_t y = (idx / dims[0]) % dims[1];
        size_t z = idx / (dims[0]*dims[1]);
        size_t sample_x = x / (1/std::cbrt(sample_rate));
        size_t sample_y = y / (1/std::cbrt(sample_rate));
        size_t sample_z = z / (1/std::cbrt(sample_rate));
        return sample_x + sample_dims[0] * sample_y + sample_dims[1] * sample_z;
    }
}

__global__ void estimate_SZp_ssim_kernel_f32(const float* oriData, size_t n_dim, size_t dims[3], float* sampleData, float* decData, 
                                             const float eb, const size_t nbEle, const size_t sample_stride)
{
    const int tid = threadIdx.x;
    const int idx = blockIdx.x * blockDim.x + tid;
    const float recipPrecision = 0.5f / eb;
    const int block_num = cmp_chunk_f32/32;
    const int start_idx = idx * cmp_chunk_f32;
    const int start_block_idx = start_idx/32;
    size_t nbBlocks = nbEle / 32;
    size_t nbSampleBlocks = (nbBlocks / sample_stride);
    size_t nbSampleEle = nbSampleBlocks * 32;
    int currQuant, lorenQuant, prevQuant, maxQuant;
    int temp_start_idx, temp_end_idx;
    int quant_chunk_idx;
    int block_idx;

    size_t sample_dims[3];
    size_t num_sample_points = 1;
    float sample_rate = 1/sample_stride;
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

    for (int j=0; j<block_num; j+=sample_stride) {
        temp_start_idx = start_idx + j*32;
        temp_end_idx = temp_start_idx + 32;
        block_idx = start_block_idx+j;
        for(int i=temp_start_idx; i<temp_end_idx; i++){
            // compute lorenzo quantization index
            float curr = oriData[i];
            currQuant = quantization_f32(curr, recipPrecision);
			float pointEstimate = currQuant * 2 * eb;

            size_t sample_idx = i / (32 * sample_stride) * 32 + (i % 32);
            // size_t sample_idx = multidim_sample_idx(i, n_dim, dims, sample_dims, sample_stride);
            if (sample_idx < nbSampleEle) {
                sampleData[sample_idx] = curr;
                decData[sample_idx] = pointEstimate;
            }

        }
    }
}