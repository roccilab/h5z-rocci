/**
 *  @file estimateCR.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <SZx/szx_float.h>
#include <SZx/sampling.h>

float szx_estimate_cr_float(const float* data, float absErrorBound, int samplingStride, size_t nbEle) {
    int blockSize = 128;
    float* radiusArray = NULL;
    float* mediusArray = NULL;
    float* buffer = NULL;
    float* data_ = (float*) malloc(sizeof(float)*nbEle);
    memcpy(data_, data, sizeof(float)*(nbEle));
    float approximateValueRange = computeRadiusBuffer_float(data_, nbEle, samplingStride, blockSize, &radiusArray, &mediusArray, &buffer);    
    size_t sumReqNbB = 0, sum_actual_leadNum = 0;
    float CR = estimateCRbasedonErrorBound_buffered_float(absErrorBound, buffer, mediusArray, radiusArray, samplingStride, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);	

    free(radiusArray);
    free(mediusArray);
    free(buffer);
    free(data_);

    return CR;
}

float szx_estimate_ssim_float(const float* data, float absErrorBound, int samplingStride, int blockSize, size_t dims[3], size_t nbEle) {
    float* data_ = (float*)malloc(sizeof(float)*nbEle);
    memcpy(data_, data, sizeof(float)*(nbEle));
    float* sample_data;
    size_t sample_dims[3];
    size_t sampleNbEle;
    size_t n_dims = 3;
    if(dims[2] == 0){
		n_dims = 2;
	}
	else if (dims[1] == 0){
	    n_dims = 1;
	}

    if(n_dims == 1){
        sample_1d_float(data_, blockSize, samplingStride, dims, &sample_data, &sampleNbEle, sample_dims);
    }
    else if (n_dims == 2){
        sample_2d_float(data_, blockSize, samplingStride, dims, &sample_data, &sampleNbEle, sample_dims);
    }
    else{ // n_dims == 3
        sample_3d_float(data_, blockSize, samplingStride, dims, &sample_data, &sampleNbEle, sample_dims);
    }

    int szx_blockSize = 128;
    size_t szx_nbBlocks = sampleNbEle / szx_blockSize;
    unsigned char* stateArray = (unsigned char*)malloc(sizeof(unsigned char)*(szx_nbBlocks+1));
    float* radiusArray = (float*)malloc(sizeof(float)*(szx_nbBlocks+1));
    float* mediusArray = (float*)malloc(sizeof(float)*(szx_nbBlocks+1));
    float approximateValueRange = computeStateMedianRadius_float(sample_data, sampleNbEle, absErrorBound, szx_blockSize, stateArray, mediusArray, radiusArray);    
    size_t sumReqNbB = 0, sum_actual_leadNum = 0;
    float SSIM = estimateSSIMbasedonErrorBound_buffered_float(absErrorBound, sample_data, mediusArray, radiusArray, szx_blockSize, sampleNbEle, &sumReqNbB, &sum_actual_leadNum, sample_dims, n_dims);	

    free(radiusArray);
    free(mediusArray);
    free(data_);
    free(sample_data);

    return SSIM;
}

float szx_estimate_psnr_float(const float* data, float absErrorBound, int samplingStride, size_t nbEle) {
    float* data_ = (float*)malloc(sizeof(float)*nbEle);
    memcpy(data_, data, sizeof(float)*(nbEle));
    float max = data[0], min = data[0];
    for(int j = 0; j < nbEle; j++){
        if(data[j] < min){
            min = data[j];
        }
        if(data[j] > max){
            max = data[j];
        }
    }

    float val_range = fabs(max-min);
    int blockSize = 128;
    float* radiusArray = NULL;
    float* mediusArray = NULL;
    float* buffer = NULL;
    float approximateValueRange = computeRadiusBuffer_float(data_, nbEle, samplingStride, blockSize, &radiusArray, &mediusArray, &buffer);    
    size_t sumReqNbB = 0, sum_actual_leadNum = 0;
    float psnr = estimatePSNRbasedonErrorBound_buffered_float(absErrorBound, buffer, val_range, mediusArray, radiusArray, samplingStride, blockSize, nbEle, &sumReqNbB, &sum_actual_leadNum);	

    free(radiusArray);
    free(mediusArray);
    free(buffer);
    free(data_);

    return psnr;
}
