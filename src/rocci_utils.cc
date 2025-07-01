#include <rocci_utils.h>
#include <string.h>

size_t computeNbEle(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1){
	size_t n = 0;
    if (r2 == 0) {
        n = r1;
    } else if (r3 == 0) {
        n = r1 * r2;
    } else if (r4 == 0) {
        n = r1 * r2 * r3;
    } else if (r5 == 0) {
        n = r1 * r2 * r3 * r4;
    } else {
        n = r1 * r2 * r3 * r4 * r5;
    }

	return n;
}

char* append_string(const char* original, const char* append) {
    size_t len1 = strlen(original);
    size_t len2 = strlen(append);
    
    char* result = (char*) malloc(len1 + len2 + 1);
    
    if (result == NULL) {
        return NULL;
    }
    
    strcpy(result, original);
    strcat(result, append);
    
    return result;
}


enum pressio_dtype rocci_dtype_to_pressio_dtype(int dataType){
	switch(dataType){
		case ROCCI_FLOAT:
			return pressio_float_dtype;
		case ROCCI_DOUBLE:
			return pressio_double_dtype;
		case ROCCI_UINT8:
			return pressio_uint8_dtype;
		case ROCCI_INT8:
			return pressio_int8_dtype;
		case ROCCI_UINT16:
			return pressio_uint16_dtype;
		case ROCCI_INT16:
			return pressio_int16_dtype;
		case ROCCI_UINT32:
			return pressio_uint32_dtype;
		case ROCCI_INT32:
			return pressio_int32_dtype;
		case ROCCI_UINT64:
			return pressio_uint64_dtype;
		case ROCCI_INT64:
			return pressio_int64_dtype;
		default:
			printf("Unrecognized data type\n");
			exit(0);
	}
}

struct pressio_compressor* get_compressor(ROCCI_Setting rocci_config, struct pressio* library){
	char* compressor_id;
	switch(rocci_config.compressorID){
		case ROCCI_SZ2:
			compressor_id = "sz";
			break;
		case ROCCI_SZ3:
			compressor_id = "sz3";
			break;
		case ROCCI_QOZ:
			compressor_id = "qoz";
			break;
		case ROCCI_ZFP:
			compressor_id = "zfp";
			break;
		case ROCCI_DR:
			compressor_id = "dr";
			break;
		case ROCCI_BG:
			compressor_id = "bg";
			break;
		case ROCCI_SZX:
			compressor_id = "szx";
			break;
		default:
			printf("Undefined compressor ID in ROCCI Configuration\n");
			exit(1);
	}

	printf("compressor: %s\n", compressor_id);
	struct pressio_compressor* comp = pressio_get_compressor(library, compressor_id);
	if(comp == NULL) {
		printf("Error: Chosen compressor, %s, is not available in libpressio. This is usually because it was not installed.\n", compressor_id);
		exit(1);
	}
	return comp;
}

// get SECRE surrogate plugin for compressor
char* get_surrogate(ROCCI_Setting rocci_config){
	char* compressor_id;
	switch(rocci_config.compressorID){
		case ROCCI_SZ3:
			compressor_id = "sz3_surrogate";
			break;
		case ROCCI_ZFP:
			compressor_id = "zfp_surrogate";
			break;
		case ROCCI_SZX:
			compressor_id = "szx_surrogate";
			break;
		case ROCCI_SZP:
			compressor_id = "szp_surrogate";
			break;
		case ROCCI_CUSZP:
			compressor_id = "cuszp_surrogate";
			break;
		case ROCCI_SPERR:
			compressor_id = "sperr_surrogate";
			break;
		default:
			return NULL;
	}

	printf("compressor surrogate: %s\n", compressor_id);
	return compressor_id;
}

struct pressio_options* get_compressor_options(ROCCI_Setting rocci_config)
{
	struct pressio_options* comp_options = pressio_options_new();
	switch(rocci_config.compressionMode)
	{
		case ABS:
			pressio_options_set_double(comp_options, "pressio:abs", rocci_config.absErrorBound);
			break;
		case REL:
			pressio_options_set_double(comp_options, "pressio:rel", rocci_config.relErrorBound);
			break;
		default:
			printf("Error: unknown compression mode\n");
			exit(0);
	}

	return comp_options;
}

struct pressio_data* create_pressio_data(int dataType, void* data, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t nDims;
	size_t nbEle;
	size_t dims[5];
	if(r2 == 0){
		nDims = 1;
		nbEle = r1;
		dims[0] = r1;
	} else if (r3 == 0){
		nDims = 2;
		nbEle = r1*r2;
		dims[0] = r2;
		dims[1] = r1;
	} else if (r4 == 0){
		nDims = 3;
		nbEle = r1*r2*r3;
		dims[0] = r3;
		dims[1] = r2;
		dims[2] = r1;
	} else if (r5 == 0){
		nDims = 4;
		nbEle = r1*r2*r3*r4;
		dims[0] = r4;
		dims[1] = r3;
		dims[2] = r2;
		dims[3] = r1;
	} else {
		nDims = 5;
		nbEle = r1*r2*r3*r4*r5;
		dims[0] = r5;
		dims[1] = r4;
		dims[2] = r3;
		dims[3] = r2;
		dims[4] = r1;
	}

	enum pressio_dtype data_type = rocci_dtype_to_pressio_dtype(dataType);
	struct pressio_data* data_buffer;
	if(data == NULL){
		float* x = (float*)malloc(sizeof(float)*nbEle);
		data_buffer = pressio_data_new_move(data_type, x, nDims, dims, NULL, NULL);
	}
	else {
		data_buffer = pressio_data_new_move(data_type, data, nDims, dims, pressio_data_libc_free_fn, NULL);
	}
	
	return data_buffer;
}