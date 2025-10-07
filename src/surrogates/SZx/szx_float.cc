#include <qcat_ssim.h>
#include <SZx/szx_float.h>
#include <math.h>

short getPrecisionReqLength_double(double precision)
{
	ldouble lbuf;
	lbuf.value = precision;
	long lvalue = lbuf.lvalue;
	
	int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
	expValue -= 1023;
//	unsigned char the1stManBit = (unsigned char)((lvalue & 0x0008000000000000) >> 51);
//	if(the1stManBit==1)
//		expValue--;
	return (short)expValue;
}

void computeReqLength_float(double realPrecision, short radExpo, int *reqLength, float *medianValue) {
    short reqExpo = getPrecisionReqLength_double(realPrecision);
    *reqLength = 9 + radExpo - reqExpo + 1; //radExpo-reqExpo == reqMantiLength
    if (*reqLength < 9)
        *reqLength = 9;
    if (*reqLength > 32) {
        *reqLength = 32;
        *medianValue = 0;
    }
}

short getExponent_float(float value)
{
	//int ivalue = floatToBigEndianInt(value);

	lfloat lbuf;
	lbuf.value = value;
	int ivalue = lbuf.ivalue;
	
	int expValue = (ivalue & 0x7F800000) >> 23;
	expValue -= 127;
	return (short)expValue;
}

// return float value with fewer bits of precision in the mantissa
float lose_mantissa_bits(float value, int bits_to_lose) {
    // Interpret the float as an integer
    unsigned int* float_bits = (unsigned int*)&value;

    // Define a mask to clear the least significant bits of the mantissa
    unsigned int mantissa_mask = bits_to_lose > 0 ? (0xFFFFFFFF >> (23 - bits_to_lose)) : 0;

    // Apply the mask to clear the specified number of bits from the mantissa
    *float_bits &= ~mantissa_mask;

    return value;
}

float computeRadiusBuffer_float(float *oriData, size_t nbEle, int samplingRate, int blockSize, float** radiusArray, float** mediusArray, float** buffer)
{
	size_t i = 0, j = 0;
	size_t offset = 0;
	size_t stride = samplingRate*blockSize;
	size_t nbBlocks = nbEle/stride;
	*radiusArray = (float*)malloc(sizeof(float)*(nbBlocks+1));
	*mediusArray = (float*)malloc(sizeof(float)*(nbBlocks+1));
	*buffer = (float*)malloc(sizeof(float)*(nbBlocks+1)*blockSize);
	float* p = *buffer;
	float g_min = oriData[offset], g_max = oriData[offset];
	for(i=0;i<nbBlocks;i++)
	{
		memcpy(p, &oriData[offset], blockSize*sizeof(float));
		float min = oriData[offset];
		float max = oriData[offset];
        for (j = 1; j < blockSize; j++) {
		float v = oriData[offset + j];
		if (min > v)
			min = v;
		else if (max < v)
			max = v;
		}
		if(g_max < max)
			g_max = max;
		if(g_min > min)
			g_min = min;
        float valueRange = max - min;
        float radius = valueRange / 2;		
		(*radiusArray)[i] = radius;
		(*mediusArray)[i] = min + radius;
		offset += stride;
		p += blockSize;
	}

	return g_max - g_min;
}

size_t computeStateMedianRadius_float(float *oriData, size_t nbEle, float absErrBound, int blockSize,
                                      unsigned char *stateArray, float *medianArray, float *radiusArray) {
    size_t nbConstantBlocks = 0;
    size_t i = 0, j = 0;
    size_t nbBlocks = nbEle / blockSize;
    size_t offset = 0;

    for (i = 0; i < nbBlocks; i++) {
        float min = oriData[offset];
        float max = oriData[offset];
        for (j = 1; j < blockSize; j++) {
            float v = oriData[offset + j];
            if (min > v)
                min = v;
            else if (max < v)
                max = v;
        }
        float valueRange = max - min;
        float radius = valueRange / 2;
        float medianValue = min + radius;

        if (radius <= absErrBound) {
            stateArray[i] = 0;
            nbConstantBlocks++;
        } else
            stateArray[i] = 1;

        //stateArray[i] = radius <= absErrBound ? 0 : 1;
        medianArray[i] = medianValue;
        radiusArray[i] = radius;
        offset += blockSize;
    }

    int remainCount = nbEle % blockSize;
    if (remainCount != 0) {
        float min = oriData[offset];
        float max = oriData[offset];
        for (j = 1; j < remainCount; j++) {
            float v = oriData[offset + j];
            if (min > v)
                min = v;
            else if (max < v)
                max = v;
        }
        float valueRange = max - min;
        float radius = valueRange / 2;
        float medianValue = min + radius;
        if (radius <= absErrBound) {
            stateArray[i] = 0;
            nbConstantBlocks++;
        } else
            stateArray[i] = 1;
        medianArray[i] = medianValue;
        radiusArray[i] = radius;
    }
    return nbConstantBlocks;
}


float estimateSSIMbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int blockSize, size_t nbEle, 
size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers, size_t dims[], size_t n_dims)
{
	size_t nbSampledBlocks = nbEle/(blockSize);
	size_t nbConstantBlocks = 0;
	*sum_actual_leadNumbers = 0;
	float metadata = 13.0*blockSize/nbEle;
	float block_cost = 33.0/8;	
	size_t i = 0, j = 0;
	*sumReqNbBytes = 0;
	float* decBuffer = (float*)malloc(sizeof(float) * nbEle);
	size_t decIndex = 0;

	float ssimSum = 0;
	for(i=0;i<nbSampledBlocks;i++) 
	{
		size_t offset = i*blockSize;
		float medianValue = medianArray[i];
		float radius = radiusArray[i];
		float *data = &buffer[offset];

//		if(i==463)
//			printf("i=%zu\n", i);
		
		if (radius <= errorBound) {
			nbConstantBlocks++;

			// simulate constant block
			for(j = 0; j < blockSize; j++){
				decBuffer[decIndex++] = medianValue;
			}

		}
		else //non-constant
		{
			int reqLength;
			short radExpo = getExponent_float(radius);
			computeReqLength_float(errorBound, radExpo, &reqLength, &medianValue);
			int reqBytesLength = reqLength / 8;
			int resiBitsLength = reqLength % 8;
			int rightShiftBits = 0;	

			if (resiBitsLength != 0) {
				rightShiftBits = 8 - resiBitsLength;				
				reqBytesLength++;
			}
			
			//printf("%d\n",reqBytesLength);
			*sumReqNbBytes+=	reqBytesLength;

			// we only do bit truncation here, so the data will simply by bitshifted
			for(j=0;j<blockSize;j++)
			{
				float data_shifted = lose_mantissa_bits(data[j], rightShiftBits);
				decBuffer[decIndex++] = data_shifted;
			}

			// *sum_actual_leadNumbers += s_actual_leadNumbers;			
		}		
	}

	// printf("NUM CONSTANT BLOCKS: %i\n", nbConstantBlocks);
	
	double estimatedSSIM;
	// printf("num_dims %i\n", n_dims);
	if(n_dims == 1){
		estimatedSSIM = SSIM_1d_windowed_float(buffer, decBuffer, dims[0], 8, 8);
	}
	else if (n_dims == 2){
		estimatedSSIM = SSIM_2d_windowed_float(buffer, decBuffer, dims[1], dims[0], 8,8, 8,8);
	}
	else { // n_dims == 3
		estimatedSSIM = SSIM_3d_windowed_float(buffer, decBuffer, dims[2], dims[1], dims[0], 8,8,8, 8,8,8);
	}

	//printf("----> sum_actual_leadNumbers=%zu, nbSampledBlocks = %zu, nbConstantBlocks = %zu, sumReqNbBytes = %zu, avgReqNbBytes=%f, avg_actual_lead=%f, p_lambda=%f\n", *sum_actual_leadNumbers, nbSampledBlocks, nbConstantBlocks, *sumReqNbBytes, avgReqNbBytes, avg_actual_lead, p_lambda);
	return estimatedSSIM;	
}

float estimatePSNRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float value_range, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, 
size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers)
{
	size_t nbSampledBlocks = nbEle/(blockSize*samplingRate); //ignored the last remainder block
	size_t nbConstantBlocks = 0;
	size_t nbNonConstantBlocks = 0;
	*sum_actual_leadNumbers = 0;
	float metadata = 13.0*blockSize/nbEle;
	float block_cost = 33.0/8;	
	size_t i = 0, j = 0;
	*sumReqNbBytes = 0;

	float maxVal = buffer[0];
	float minVal = buffer[0];
	float sq_err = 0;

	for(i=0;i<nbSampledBlocks;i++) //ignored the last remainder block
	{
		size_t offset = i*blockSize;
		float medianValue = medianArray[i];
		float radius = radiusArray[i];

		float *data = &buffer[offset];

		// if(minVal > min){
		//     minVal = min;
		// }
		// if(maxVal < max){
		//     maxVal = max;
		// }

//		if(i==463)
//			printf("i=%zu\n", i);
		
		if (radius <= errorBound) { // constant block
			nbConstantBlocks++;

			// update mse
			for(int k = 0; k < blockSize; k++){
				sq_err += pow((data[k] - medianValue), 2);
			}

		}
		else //non-constant
		{
			nbNonConstantBlocks++;
			int reqLength;
			short radExpo = getExponent_float(radius);
			computeReqLength_float(errorBound, radExpo, &reqLength, &medianValue);
			int reqBytesLength = reqLength / 8;
			int resiBitsLength = reqLength % 8;
			int rightShiftBits = 0;	

			if (resiBitsLength != 0) {
				rightShiftBits = 8 - resiBitsLength;				
				reqBytesLength++;
			}
			
			//printf("%d\n",reqBytesLength);
			*sumReqNbBytes+=	reqBytesLength;
			
			register lfloat lfBuf_pre;
			register lfloat lfBuf_cur;
			lfBuf_pre.ivalue = 0;
			register unsigned char leadingNum = 0;			
			//int leadingNum_Array[128];
			int s_actual_leadNumbers = 0;

			// we only do bit truncation here, so the data will simply by bitshifted
			for(j=0;j<blockSize;j++)
			{
				lfBuf_cur.value = data[j] - medianValue;

				lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

				lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;
				
				leadingNum = 0;
				if (lfBuf_pre.ivalue >> 8 == 0)
					leadingNum = 3;
				else if (lfBuf_pre.ivalue >> 16 == 0)
					leadingNum = 2;
				else if (lfBuf_pre.ivalue >> 24 == 0)
					leadingNum = 1;
					
				// leadingNum_Array[j] = leadingNum;    
				if(leadingNum >= reqBytesLength)
					s_actual_leadNumbers += reqBytesLength;
				else
					s_actual_leadNumbers += leadingNum;     
					
				lfBuf_pre = lfBuf_cur;	            

				float data_shifted = lose_mantissa_bits(data[j], rightShiftBits);
				sq_err += pow(errorBound, 2); //pow((data_shifted - data[j]), 2);
			}

			*sum_actual_leadNumbers += s_actual_leadNumbers;			
		}		
	}
	
	float avgReqNbBytes = 1.0*(*sumReqNbBytes)/(nbSampledBlocks - nbConstantBlocks);
	float avg_actual_lead = 1.0*(*sum_actual_leadNumbers)/(nbSampledBlocks - nbConstantBlocks);
	
	float p_lambda = 1.0*nbConstantBlocks/nbSampledBlocks;
	
	size_t samplenbEle = blockSize * nbSampledBlocks;
	float mse = sq_err / nbEle;//(samplenbEle);
	// printf("SQERR %f, MSE %.10f, MAX %f, SampleNBELE %i\n", sq_err, mse, maxVal, samplenbEle);
	// for the case where mse == 0 we add a small constant to avoid numerical errors
	float eps = 1e-16;
	// float value_range = maxVal - minVal;
	float estimatedPSNR = -20.0*log10((sqrt(mse) / value_range) + eps);

	//printf("----> sum_actual_leadNumbers=%zu, nbSampledBlocks = %zu, nbConstantBlocks = %zu, sumReqNbBytes = %zu, avgReqNbBytes=%f, avg_actual_lead=%f, p_lambda=%f\n", *sum_actual_leadNumbers, nbSampledBlocks, nbConstantBlocks, *sumReqNbBytes, avgReqNbBytes, avg_actual_lead, p_lambda);
	return estimatedPSNR;	
}

float estimateCRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, 
size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers)
{
	size_t nbSampledBlocks = nbEle/(blockSize*samplingRate); //ignored the last remainder block
	size_t nbConstantBlocks = 0;
	*sum_actual_leadNumbers = 0;
	float metadata = 13.0*blockSize/nbEle;
	float block_cost = 33.0/8;	
	size_t i = 0, j = 0;
	*sumReqNbBytes = 0;
	for(i=0;i<nbSampledBlocks;i++) //ignored the last remainder block
	{
		size_t offset = i*blockSize;
		float medianValue = medianArray[i];
		float radius = radiusArray[i];
		float *data = &buffer[offset];
//		if(i==463)
//			printf("i=%zu\n", i);
		
		if (radius <= errorBound) {
			nbConstantBlocks++;
		}
		else //non-constant
		{
			int reqLength;
			short radExpo = getExponent_float(radius);
			computeReqLength_float(errorBound, radExpo, &reqLength, &medianValue);
			int reqBytesLength = reqLength / 8;
			int resiBitsLength = reqLength % 8;
			int rightShiftBits = 0;	

			if (resiBitsLength != 0) {
				rightShiftBits = 8 - resiBitsLength;				
				reqBytesLength++;
			}
			
			//printf("%d\n",reqBytesLength);
			*sumReqNbBytes+=	reqBytesLength;
			
			register lfloat lfBuf_pre;
			register lfloat lfBuf_cur;
			lfBuf_pre.ivalue = 0;
			register unsigned char leadingNum = 0;			
			//int leadingNum_Array[128];
			int s_actual_leadNumbers = 0;
			for(j=0;j<blockSize;j++)
			{
				lfBuf_cur.value = data[j] - medianValue;

				lfBuf_cur.ivalue = lfBuf_cur.ivalue >> rightShiftBits;

				lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;
				
				leadingNum = 0;
				if (lfBuf_pre.ivalue >> 8 == 0)
					leadingNum = 3;
				else if (lfBuf_pre.ivalue >> 16 == 0)
					leadingNum = 2;
				else if (lfBuf_pre.ivalue >> 24 == 0)
					leadingNum = 1;
					
				// leadingNum_Array[j] = leadingNum;    
				if(leadingNum >= reqBytesLength)
					s_actual_leadNumbers += reqBytesLength;
				else
					s_actual_leadNumbers += leadingNum;     
					
				lfBuf_pre = lfBuf_cur;	               
			}
			
			*sum_actual_leadNumbers += s_actual_leadNumbers;			
		}		
	}
	
	float avgReqNbBytes = 1.0*(*sumReqNbBytes)/(nbSampledBlocks - nbConstantBlocks);
	float avg_actual_lead = 1.0*(*sum_actual_leadNumbers)/(nbSampledBlocks - nbConstantBlocks);
	
	float p_lambda = 1.0*nbConstantBlocks/nbSampledBlocks;

	if(nbSampledBlocks == nbConstantBlocks) {
		// avoid Nan
		avgReqNbBytes = sizeof(float);
		avg_actual_lead = 0.0;
	}
	
	float estimatedCR = 4*blockSize/(metadata + block_cost+(1 + (0.25+avgReqNbBytes)*blockSize - avg_actual_lead)*(1 - p_lambda));
	//printf("----> sum_actual_leadNumbers=%zu, nbSampledBlocks = %zu, nbConstantBlocks = %zu, sumReqNbBytes = %zu, avgReqNbBytes=%f, avg_actual_lead=%f, p_lambda=%f\n", *sum_actual_leadNumbers, nbSampledBlocks, nbConstantBlocks, *sumReqNbBytes, avgReqNbBytes, avg_actual_lead, p_lambda);
	return estimatedCR;	
}

