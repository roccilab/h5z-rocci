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
#include <rocci.h>
#include <rocci_utils.h>

int versionNumber[4] = {ROCCI_VER_MAJOR,ROCCI_VER_MINOR,ROCCI_VER_BUILD,ROCCI_VER_REVISION};

const char* COMPRESSOR_STR[] = {"SZ2", "SZ3", "QOZ", "ZFP", "DR", "BG", "SZX"};
int COMPRESSOR_OPTIONS[] = {ROCCI_SZ2, ROCCI_SZ3, ROCCI_QOZ, ROCCI_ZFP, ROCCI_DR, ROCCI_BG, ROCCI_SZX};

const char* CMP_MODE_STR[] = {"ABS", "REL", "VR_REL", "ABS_AND_REL", "ABS_OR_REL", "PSNR", "NORM", "FIX_RATE"};
int CMP_MODE_OPTIONS[] = {ABS, REL, VR_REL, ABS_AND_REL, ABS_OR_REL, PSNR, NORM, FIX_RATE};

ROCCI_Setting rocci_config;

static int rocci_ini_config_handler(void* user, const char* section, const char* name,
                   const char* value)
{
    ROCCI_Setting* config = (ROCCI_Setting*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

	if(MATCH("GlobalSettings", "CompressorID")){
		int nCompressors = sizeof(COMPRESSOR_STR) / sizeof(char*);
		for(int i = 0; i < nCompressors; i++){
			if(strcmp(value, COMPRESSOR_STR[i]) == 0){
				config->compressorID = COMPRESSOR_OPTIONS[i];
			}
		}
	}
	else if(MATCH("GlobalSettings", "CompressionMode")){
		int nModes = sizeof(CMP_MODE_STR) / sizeof(char*);
		for(int i = 0; i < nModes; i++){
			if(strcmp(value, CMP_MODE_STR[i]) == 0){
				config->compressionMode = CMP_MODE_OPTIONS[i];
			}
		}
	}
	else if(MATCH("GlobalSettings", "AbsErrorBound")) {
		config->absErrorBound = atof(value);
	}
	else if(MATCH("GlobalSettings", "RelErrorBound")) {
		config->relErrorBound = atof(value);
	}
	else if(MATCH("GlobalSettings", "Precision")) {
		config->precision = atoi(value);
	}
	else if(MATCH("GlobalSettings", "FixedCompressionRate")) {
		config->rate = atof(value);
	} 
	else {
		printf("Unknown section or name: %s/%s\n", section, name);
        return 0;
    }
    return 1;
}

void ROCCI_Init(char* cfgFile)
{
	if(cfgFile == NULL){
		rocci_config.compressorID = ROCCI_SZ3;
		rocci_config.compressionMode = ABS;
		rocci_config.absErrorBound = 0.01;
		rocci_config.relErrorBound = 0;
		rocci_config.precision = 0;
		rocci_config.rate = 0;
		printf("Warning: No config file found, using defaults\n");
		return;
	}

	printf("Reading %s...\n", cfgFile);
	if (ini_parse(cfgFile, rocci_ini_config_handler, &rocci_config) < 0) {
        printf("Unable to load config file: %s\n", cfgFile);
        exit(0);
    }

	printf("Succesfully loaded config: %s\n", cfgFile);
	// printf("CONFIG: %i %i %f %f %i %f\n", rocci_config.compressorID, rocci_config.compressionMode, rocci_config.absErrorBound, rocci_config.relErrorBound, rocci_config.precision, rocci_config.rate);
}

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
	struct pressio* library = pressio_instance();
	struct pressio_compressor* comp = get_compressor(rocci_config, library);
	struct pressio_data* input_data = create_pressio_data(dataType, data, r5, r4, r3, r2, r1);
    struct pressio_data* compressed = pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
	struct pressio_options* comp_options = get_compressor_options(rocci_config);
	if(pressio_compressor_set_options(comp, comp_options)) {
		fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
		return NULL;
	}
  	pressio_options_free(comp_options);

	if(pressio_compressor_compress(comp, input_data, compressed)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
		return NULL;
    }

	void* output_data = pressio_data_ptr(compressed, outSize);
	pressio_compressor_release(comp);
	pressio_release(library);

	return (unsigned char*)output_data;
}

unsigned char* ROCCI_compress(int dataType, void* data, size_t* outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	struct pressio* library = pressio_instance();
	struct pressio_compressor* comp = get_compressor(rocci_config, library);
	struct pressio_data* input_data = create_pressio_data(dataType, data, r5, r4, r3, r2, r1);
    struct pressio_data* compressed = pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
	struct pressio_options* comp_options = get_compressor_options(rocci_config);
	if(pressio_compressor_set_options(comp, comp_options)) {
		fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
		return NULL;
	}
  	pressio_options_free(comp_options);

	if(pressio_compressor_compress(comp, input_data, compressed)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
		return NULL;
    }

	void* output_data = pressio_data_ptr(compressed, outSize);
	pressio_compressor_release(comp);
	pressio_release(library);

	return (unsigned char*)output_data;
}

void* ROCCI_decompress(int dataType, char *buf, size_t nbytes, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t nbEle = computeNbEle(r5, r4, r3, r2, r1);

	struct pressio* library = pressio_instance();
  	struct pressio_compressor* comp = get_compressor(rocci_config, library);
	size_t dim[1] = {nbytes};
    struct pressio_data* compressed_data = pressio_data_new_move(pressio_byte_dtype, buf, 1, dim, NULL, NULL);
	struct pressio_data* output = create_pressio_data(dataType, NULL, r5, r4, r3, r2, r1);

	if(pressio_compressor_decompress(comp, compressed_data, output)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
	  exit(0);
    }

	size_t outSize = 0;
	float* output_data = pressio_data_ptr(output, &outSize);
	pressio_compressor_release(comp);
	pressio_release(library);

	return output_data;
}


