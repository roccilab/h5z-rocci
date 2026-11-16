#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <SZx/szx_float.h>
#include <SZx/sampling.h>

// Inline non-buffered CR estimation (from SZx source, exit(0) removed).
// No malloc/memcpy/precomputation overhead.
static float estimate_cr_direct(float errorBound, const float* data, int blockSize, size_t nbEle,
                                int samplingRate) {
    float metadata = 13.0f * blockSize / nbEle;
    float block_cost = 33.0f / 8;
    size_t nbBlocks = nbEle / blockSize;
    size_t nbConstantBlocks = 0;
    size_t sum_actual_leadNumbers = 0;
    if (samplingRate < 1) samplingRate = 1;
    size_t nbBlocks_sampled = nbBlocks / samplingRate;
    size_t sumReqNbBytes = 0;

    for (size_t i = 0; i < nbBlocks; i += samplingRate) {
        size_t offset = i * blockSize;
        float min = data[offset], max = data[offset];
        for (size_t j = 1; j < (size_t)blockSize; j++) {
            float v = data[offset + j];
            if (min > v) min = v;
            if (max < v) max = v;
        }
        float radius = (max - min) / 2;

        if (radius <= errorBound) {
            nbConstantBlocks++;
        } else {
            short radExpo = getExponent_float(radius);
            float medianValue = min + radius;
            int reqLength;
            computeReqLength_float(errorBound, radExpo, &reqLength, &medianValue);
            int reqBytesLength = reqLength / 8;
            int resiBitsLength = reqLength % 8;
            int rightShiftBits = 0;
            if (resiBitsLength != 0) {
                rightShiftBits = 8 - resiBitsLength;
                reqBytesLength++;
            }
            sumReqNbBytes += reqBytesLength;

            typedef union { float value; unsigned int ivalue; } lfloat;
            lfloat lfBuf_pre, lfBuf_cur;
            lfBuf_pre.ivalue = 0;
            int s_actual_leadNumbers = 0;
            for (size_t j = 0; j < (size_t)blockSize; j++) {
                lfBuf_cur.value = data[offset + j] - medianValue;
                lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;
                lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;
                unsigned char leadingNum = 0;
                if (lfBuf_pre.ivalue >> 8 == 0) leadingNum = 3;
                else if (lfBuf_pre.ivalue >> 16 == 0) leadingNum = 2;
                else if (lfBuf_pre.ivalue >> 24 == 0) leadingNum = 1;
                if (leadingNum >= reqBytesLength)
                    s_actual_leadNumbers += reqBytesLength;
                else
                    s_actual_leadNumbers += leadingNum;
                lfBuf_pre = lfBuf_cur;
            }
            sum_actual_leadNumbers += s_actual_leadNumbers;
        }
    }

    float avgReqNbBytes = (nbBlocks_sampled > nbConstantBlocks)
        ? 1.0f * sumReqNbBytes / (nbBlocks_sampled - nbConstantBlocks) : 4.0f;
    float avg_actual_lead = (nbBlocks_sampled > nbConstantBlocks)
        ? 1.0f * sum_actual_leadNumbers / (nbBlocks_sampled - nbConstantBlocks) : 0.0f;
    float p_lambda = 1.0f * nbConstantBlocks / (nbBlocks_sampled > 0 ? nbBlocks_sampled : 1);

    return 4.0f * blockSize / (metadata + block_cost +
        (1 + (0.25f + avgReqNbBytes) * blockSize - avg_actual_lead) * (1 - p_lambda));
}

float szx_estimate_cr_float(const float* data, float absErrorBound, int samplingStride, size_t nbEle) {
    return estimate_cr_direct(absErrorBound, data, 128, nbEle, samplingStride);
}

float szx_estimate_ssim_float(const float* data, float absErrorBound, int samplingStride, int blockSize, size_t dims[3], size_t nbEle) {
    float* data_ = (float*)malloc(sizeof(float)*nbEle);
    memcpy(data_, data, sizeof(float)*(nbEle));
    float* sample_data;
    size_t sample_dims[3];
    size_t sampleNbEle;
    size_t n_dims = 3;
    if(dims[2] == 0) n_dims = 2;
    else if (dims[1] == 0) n_dims = 1;

    // SSIM sampling block size must not exceed any dimension (else that axis
    // gets 0 sampled blocks -> sample_dims=0 -> SSIM window fails).
    // Cap at smallest dim; keep >= 8 for the 8x8x8 SSIM window.
    size_t min_dim = dims[0];
    for (int a = 1; a < (int)n_dims; a++)
        if (dims[a] < min_dim) min_dim = dims[a];
    int ssim_bs = blockSize;
    if ((size_t)ssim_bs > min_dim) ssim_bs = (int)min_dim;
    if (ssim_bs < 8) ssim_bs = 8;

    if(n_dims == 1)
        sample_1d_float(data_, ssim_bs, samplingStride, dims, &sample_data, &sampleNbEle, sample_dims);
    else if (n_dims == 2)
        sample_2d_float(data_, ssim_bs, samplingStride, dims, &sample_data, &sampleNbEle, sample_dims);
    else
        sample_3d_float(data_, ssim_bs, samplingStride, dims, &sample_data, &sampleNbEle, sample_dims);

    int szx_blockSize = 128;
    size_t szx_nbBlocks = sampleNbEle / szx_blockSize;
    unsigned char* stateArray = (unsigned char*)malloc(sizeof(unsigned char)*(szx_nbBlocks+1));
    float* radiusArray = (float*)malloc(sizeof(float)*(szx_nbBlocks+1));
    float* mediusArray = (float*)malloc(sizeof(float)*(szx_nbBlocks+1));
    computeStateMedianRadius_float(sample_data, sampleNbEle, absErrorBound, szx_blockSize, stateArray, mediusArray, radiusArray);
    size_t sumReqNbB = 0, sum_actual_leadNum = 0;
    float SSIM = estimateSSIMbasedonErrorBound_buffered_float(absErrorBound, sample_data, mediusArray, radiusArray, szx_blockSize, sampleNbEle, &sumReqNbB, &sum_actual_leadNum, sample_dims, n_dims);

    free(radiusArray); free(mediusArray); free(stateArray); free(data_); free(sample_data);
    return SSIM;
}

float szx_estimate_psnr_float(const float* data, float absErrorBound, int samplingStride, size_t nbEle) {
    float max = data[0], min = data[0];
    for(size_t j = 1; j < nbEle; j++){
        if(data[j] < min) min = data[j];
        if(data[j] > max) max = data[j];
    }
    float val_range = fabs(max-min);

    int blockSize = 128;
    float* radiusArray = NULL, *mediusArray = NULL, *buffer = NULL;
    float* data_ = (float*) malloc(sizeof(float)*nbEle);
    memcpy(data_, data, sizeof(float)*(nbEle));
    computeRadiusBuffer_float(data_, nbEle, samplingStride, blockSize, &radiusArray, &mediusArray, &buffer);
    size_t sumReqNbB = 0, sum_actual_leadNum = 0;
    float psnr = estimatePSNRbasedonErrorBound_buffered_float(absErrorBound, buffer, val_range, mediusArray, radiusArray, samplingStride, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);

    free(radiusArray); free(mediusArray); free(buffer); free(data_);
    return psnr;
}
