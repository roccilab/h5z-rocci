/**
 *  @file rocci.c
 *  @author Sheng Di
 *  @date July, 2022
 *  @brief 
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "rocci.h"

int versionNumber[4] = {ROCCI_VER_MAJOR,ROCCI_VER_MINOR,ROCCI_VER_BUILD,ROCCI_VER_REVISION};

int rocciFidelity_multiFields(float** oriData, float** decData, float** fidelity, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	return 0;
}

int rocciFidelity_singleField(float* oriData, float* decData, float** fidelity, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	return 0;
}

int roccifastSearchBestSetting (int metric, int compressorID, ROCCI_Target input, ROCCI_Setting* result)
{
	return 0;
}

unsigned char* ROCCI_compress_args(int dataType, void* data, size_t* outSize, int error_mode, double abs_error, double rel_error, double pw_rel_error, double psnr, double ratio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	return NULL;
}

unsigned char* ROCCI_compress(int dataType, void* data, size_t* outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	return NULL;
}
