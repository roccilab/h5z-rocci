inline void computeReqLength_float(double realPrecision, short radExpo, int *reqLength, float *medianValue) {
    short reqExpo = getPrecisionReqLength_double(realPrecision);
    *reqLength = 9 + radExpo - reqExpo + 1; //radExpo-reqExpo == reqMantiLength
    if (*reqLength < 9)
        *reqLength = 9;
    if (*reqLength > 32) {
        *reqLength = 32;
        *medianValue = 0;
    }
}

inline short getExponent_float(float value)
{
	//int ivalue = floatToBigEndianInt(value);

	lfloat lbuf;
	lbuf.value = value;
	int ivalue = lbuf.ivalue;
	
	int expValue = (ivalue & 0x7F800000) >> 23;
	expValue -= 127;
	return (short)expValue;
}

float lose_mantissa_bits(float value, int bits_to_lose);
float computeRadiusBuffer_float(float *oriData, size_t nbEle, int samplingRate, int blockSize, float** radiusArray, float** mediusArray, float** buffer);
size_t computeStateMedianRadius_float(float *oriData, size_t nbEle, float absErrBound, int blockSize, unsigned char *stateArray, float *medianArray, float *radiusArray);
float estimateSSIMbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int blockSize, size_t nbEle, size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers, size_t dims[], size_t n_dims);
float estimatePSNRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float value_range, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers);
float estimateCRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers);

