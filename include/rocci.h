/**
 *  @file rocci.h
 *  @author Sheng Di
 *  @date April, 2022
 *  @brief Header file for the whole compressor.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ROCCI_H
#define _ROCCI_H

#include <stdio.h>
#include <stdint.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>      /* For gettimeofday(), in microseconds */
#endif
#include <time.h>          /* For time(), in seconds */

#include <rocci_defines.h>
//#include <libpressio.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ROCCI_SZ2 0
#define ROCCI_SZ3 1
#define ROCCI_QOZ 2
#define ROCCI_ZFP 3
#define ROCCI_DR 4
#define ROCCI_BG 5
#define ROCCI_SZX 6

#define ROCCI_CR_METRIC 0
#define ROCCI_PSNR_METRIC 1
#define ROCCI_SSIM_METRIC 2
#define ROCCI_AC_METRIC 3

typedef struct ROCCI_Target
{
	int metric; //specify the metric
	double targetValue_lowbound;
	double targetValue_upperbound;
} ROCCI_Target;

typedef struct ROCCI_Setting
{
	int compressorID;
	int compressionMode;
	double absErrorBound;
	double relErrorBound;
	int precision; 
	double rate;
} ROCCI_Setting;

void ROCCI_Init(char* cfgFile);
int rocciFidelity_multiFields(float** oriData, float** decData, float** fidelity, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int rocciFidelity_singleField(float* oriData, float* decData, float** fidelity, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int roccifastSearchBestSetting (int metric, int compressorID, ROCCI_Target input, ROCCI_Setting* result);
//int constuctRuntimeCompressor (ROCCI_Setting* input , pressio_compressor* output);

unsigned char* ROCCI_compress_args(int dataType, void* data, size_t* outSize, int error_mode, double abs_error, double rel_error, double pw_rel_error, double psnr, double ratio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
unsigned char* ROCCI_compress(int dataType, void* data, size_t* outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void* ROCCI_decompress(int dataType, char *buf, size_t nbytes, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ROCCI_H  ----- */
