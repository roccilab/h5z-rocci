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
#include <math.h>
#include <rocci.h>
#include <rocci_utils.h>
#include <pressio_search_defines.h>

ROCCI_Setting rocci_config;
ROCCI_Target  rocci_target;


const int versionNumber[4] = {ROCCI_VER_MAJOR,ROCCI_VER_MINOR,ROCCI_VER_BUILD,ROCCI_VER_REVISION};

const char* COMPRESSOR_STR[] = {"SZ2", "SZ3", "QOZ", "ZFP", "DR", "BG", "SZX"};
const int COMPRESSOR_OPTIONS[] = {ROCCI_SZ2, ROCCI_SZ3, ROCCI_QOZ, ROCCI_ZFP, ROCCI_DR, ROCCI_BG, ROCCI_SZX};

const char* CMP_MODE_STR[] = {"ABS", "REL", "VR_REL", "ABS_AND_REL", "ABS_OR_REL", "PSNR", "NORM", "FIX_RATE"};
const int CMP_MODE_OPTIONS[] = {ABS, REL, VR_REL, ABS_AND_REL, ABS_OR_REL, PSNR, NORM, FIX_RATE};

const char* METRIC_STR[] = {"CR", "PSNR", "SSIM"};
const int METRIC_OPTIONS[] = {ROCCI_CR_METRIC, ROCCI_PSNR_METRIC, ROCCI_SSIM_METRIC};


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
	else if(strcmp("SECRESettings", section) == 0) {
		// skip
		return 1;
	}
	else {
		printf("Unknown section or name: %s/%s\n", section, name);
        return 0;
    }
    return 1;
}

static int secre_ini_config_handler(void* user, const char* section, const char* name,
                   const char* value)
{
    ROCCI_Target* config = (ROCCI_Target*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

	if(MATCH("SECRESettings", "UseSECRE")){
		config->do_opt = atoi(value);
	}
	else if(MATCH("SECRESettings", "SECREMetric")){
		int nMetrics = sizeof(METRIC_STR) / sizeof(char*);
		for(int i = 0; i < nMetrics; i++){
			if(strcmp(value, METRIC_STR[i]) == 0){
				config->metric = METRIC_OPTIONS[i];
			}
		}
	}
	else if(MATCH("SECRESettings", "SECRETarget")){
		config->targetValue = atof(value);
	}
	else if(MATCH("SECRESettings", "SECRETolerance")){
		config->tolerance = atof(value);
	}
	else if(strcmp("GlobalSettings", section) == 0) {
		// skip
		return 1;
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

		rocci_target.metric = ROCCI_CR_METRIC;
		rocci_target.targetValue = 60.0;
		rocci_target.tolerance = 10.0;
		rocci_target.do_opt = false;

		printf("Warning: No config file found, using defaults\n");
		return;
	}

	printf("Reading %s...\n", cfgFile);
	if (ini_parse(cfgFile, rocci_ini_config_handler, &rocci_config) < 0) {
        printf("Unable to load config file: %s\n", cfgFile);
        exit(0);
    }

	if (ini_parse(cfgFile, secre_ini_config_handler, &rocci_target) < 0) {
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

int roccifastSearchBestSetting(struct pressio_data* input_data, ROCCI_Target input, ROCCI_Setting* result)
{
	struct pressio* library = pressio_instance();
	char* surrogate_id = get_surrogate(*result, library);
	if(surrogate_id == NULL) {
		printf("Warning: No surrogate available for the given compressor, skipping parameter optimization.\n");
		return 0;
	}

	char* metric_str;
	switch(input.metric) {
		case ROCCI_CR_METRIC:
			metric_str = "cr";
			break;
		case ROCCI_PSNR_METRIC:
			metric_str = "psnr";
			break;
		case ROCCI_SSIM_METRIC:
			metric_str = "ssim";
			break;
		default:
			printf("Error: Invalid search metric\n");
			exit(0);
	}

	// compute value range
	float* data = (float*) pressio_data_ptr(input_data, NULL);
	size_t nbEle = pressio_data_num_elements(input_data);
	float min = data[0], max = data[0];
	for(size_t i = 1; i < nbEle; i++) {
		float curr = data[i];
		if(curr < min) min = curr;
		if(curr > max) max = curr;
	}
	float valRange = fabs(max - min);

	char* surrogate_prefix = append_string(surrogate_id, ":");
	char* surrogate_metric_pre = append_string(surrogate_prefix, "metric");
	char* surrogate_metric_str = append_string(surrogate_prefix, metric_str);

	// set up opt meta compressor
	struct pressio_compressor* comp = pressio_get_compressor(library, "opt");
	if(comp == NULL) {
		fprintf(stderr, "failed to load compressor: %s\n", pressio_error_msg(library));
		fprintf(stderr, "check to ensure that liblibpressio_opt.so was loaded into the binary using readelf\n");
		exit(1);
	}

	// set search options
	struct pressio_options* search_configuration = pressio_options_new();
	pressio_options_set_integer(search_configuration, "opt:do_decompress", 0);
	// optimize selected surrogate model
	pressio_options_set_string(search_configuration, "opt:compressor", surrogate_id);
	// set target values for metric
	pressio_options_set_uinteger(search_configuration, "opt:objective_mode", pressio_search_mode_target);
	pressio_options_set_double(search_configuration, "opt:target", input.targetValue);
	pressio_options_set_double(search_configuration, "opt:global_rel_tolerance", input.tolerance);
	// use FraZ for search
	pressio_options_set_string(search_configuration, "opt:search", "fraz");
	// set metric info for surrogate
	pressio_options_set_string(search_configuration, surrogate_metric_pre, metric_str);
	pressio_options_set_uinteger(search_configuration, "opt:max_iterations", 100);
	pressio_options_set_string(search_configuration, "opt:search_metrics", "progress_printer");

	// optimize the absolute error bound
	const char* inputs[] = {"pressio:abs"};
	double lower_bound[] = {0.0};
	double upper_bound[] = {(double) valRange};
	size_t bound_dims[] = {1};
	struct pressio_data* lower_bound_data = pressio_data_new_nonowning(pressio_double_dtype, lower_bound, 1, bound_dims);
	struct pressio_data* upper_bound_data = pressio_data_new_nonowning(pressio_double_dtype, upper_bound, 1, bound_dims);
	pressio_options_set_strings(search_configuration, "opt:inputs", sizeof(inputs)/sizeof(inputs[0]), inputs);
	pressio_options_set_data(search_configuration, "opt:lower_bound", lower_bound_data);
  	pressio_options_set_data(search_configuration, "opt:upper_bound", upper_bound_data);

	// use the surrogate model's estimated metric value for optimization
	const char* outputs[] = {surrogate_metric_str};
	pressio_options_set_strings(search_configuration, "opt:output", sizeof(outputs)/sizeof(outputs[0]), outputs);

	if(pressio_compressor_set_options(comp, search_configuration)) {
		fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
		return pressio_compressor_error_code(comp);
	}
  	pressio_options_free(search_configuration);

	// run optimization
	struct pressio_data* compressed = pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
	if(pressio_compressor_compress(comp, input_data, compressed)) {
      fprintf(stderr, "%s\n", pressio_compressor_error_msg(comp));
      return pressio_compressor_error_code(comp);
    }

	// get final parameter settings
	struct pressio_options* metrics_results = pressio_compressor_get_metrics_results(comp);
	struct pressio_data* eb = pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
	pressio_options_get_data(metrics_results, "opt:input", &eb);
	pressio_options_free(metrics_results);

	// set in param configuration
	size_t out_bytes;
	result->absErrorBound = ((double *) pressio_data_ptr(eb, &out_bytes))[0];
	result->compressionMode = ABS;

	return 0;
}

unsigned char* ROCCI_compress_args(int dataType, void* data, size_t* outSize, int error_mode, double abs_error, double rel_error, double pw_rel_error, double psnr, double ratio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	struct pressio* library = pressio_instance();
	struct pressio_compressor* comp = get_compressor(rocci_config, library);
	struct pressio_data* input_data = create_pressio_data(dataType, data, r5, r4, r3, r2, r1);
    struct pressio_data* compressed = pressio_data_new_empty(pressio_byte_dtype, 0, NULL);

	// optimize compressor settings if enabled
	if(rocci_target.do_opt) {
		if(roccifastSearchBestSetting(input_data, rocci_target, &rocci_config) > 0) {
			exit(1);
		}
	}

	// set compressor options
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

	// optimize compressor settings if enabled
	if(rocci_target.do_opt) {
		if(roccifastSearchBestSetting(input_data, rocci_target, &rocci_config) > 0) {
			exit(1);
		}
	}

	// set compressor options
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
	float* output_data = (float*)pressio_data_ptr(output, &outSize);
	pressio_compressor_release(comp);
	pressio_release(library);

	return output_data;
}


