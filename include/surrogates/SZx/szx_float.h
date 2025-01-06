#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef union lfloat {
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
    uint16_t int16[2];
} lfloat;

typedef union ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ldouble;

short getPrecisionReqLength_double(double precision);
void computeReqLength_float(double realPrecision, short radExpo, int *reqLength, float *medianValue);
short getExponent_float(float value);

float lose_mantissa_bits(float value, int bits_to_lose);
float computeRadiusBuffer_float(float *oriData, size_t nbEle, int samplingRate, int blockSize, float** radiusArray, float** mediusArray, float** buffer);
size_t computeStateMedianRadius_float(float *oriData, size_t nbEle, float absErrBound, int blockSize, unsigned char *stateArray, float *medianArray, float *radiusArray);
float estimateSSIMbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int blockSize, size_t nbEle, size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers, size_t dims[], size_t n_dims);
float estimatePSNRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float value_range, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers);
float estimateCRbasedonErrorBound_buffered_float(float errorBound, float* buffer, float* medianArray, float* radiusArray, int samplingRate, int blockSize, size_t nbEle, size_t *sumReqNbBytes, size_t *sum_actual_leadNumbers);

#ifdef __cplusplus
}
#endif