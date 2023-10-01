/**
 *  @file H5Z_ROCCI.c
 *  @author Sheng Di
 *  @date July, 2022
 *  @brief ROCCI filter for HDF5
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "H5Z_ROCCI.h"
#include "H5PLextern.h"
#include <stdint.h>

//rocci_params* conf_params = NULL;
int load_conffile_flag = 1; //set to load configuration file at the beginning by user

int init_rocci_flag = 0; //0 means 'not yet', 1 means 'already loaded'
char cfgFile[256] = "rocci.config"; 

static int h5z_rocci_was_registered = 0;

const H5Z_class2_t H5Z_ROCCI[1] = {{
	H5Z_CLASS_T_VERS,              /* H5Z_class_t version */
	H5Z_FILTER_ROCCI, /* Filter id number */
	1,              /* encoder_present flag (set to true) */
	1,              /* decoder_present flag (set to true) */
	"ROCCI compressor/decompressor for floating-point data.", /* Filter name for debugging */
	NULL,                          /* The "can apply" callback */
	H5Z_rocci_set_local,                          /* The "set local" callback */
	H5Z_filter_rocci,   /* The actual filter function */
}};

H5PL_type_t H5PLget_plugin_type(void) {return H5PL_TYPE_FILTER;}
const void *H5PLget_plugin_info(void) {return H5Z_ROCCI;}

void ROCCI_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
int* error_bound_mode, double* abs_error, double* rel_error, double* pw_rel_error, double* psnr, double* fixedRatio)
{
	ROCCI_cdArrayToMetaData(cd_nelmts, cd_values, dimSize, dataType, r5, r4, r3, r2, r1);
	int dim = *dimSize;
	int k = dim==1?4:dim+2;
	unsigned char b[8];
	int32ToBytes_bigEndian(b, cd_values[k++]);
	*error_bound_mode = bytesToInt32_bigEndian(b);
	int32ToBytes_bigEndian(b, cd_values[k++]);
	int32ToBytes_bigEndian(b+4, cd_values[k++]);
	*abs_error = bytesToDouble(b);
	int32ToBytes_bigEndian(b, cd_values[k++]);
	int32ToBytes_bigEndian(b+4, cd_values[k++]);
	*rel_error = bytesToDouble(b);	
	int32ToBytes_bigEndian(b, cd_values[k++]);
	int32ToBytes_bigEndian(b+4, cd_values[k++]);
	*pw_rel_error = bytesToDouble(b);
	int32ToBytes_bigEndian(b, cd_values[k++]);
	int32ToBytes_bigEndian(b+4, cd_values[k++]);
	*psnr = bytesToDouble(b);
	int32ToBytes_bigEndian(b, cd_values[k++]);
	int32ToBytes_bigEndian(b+4, cd_values[k++]);
	*fixedRatio = bytesToDouble(b);	
}

/**
 * to be used in decompression and compression, inside the H5Z_filter_rocci().
 * */
void ROCCI_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1)
{
	//printf("cd_nelmts=%zu\n", cd_nelmts); 
	assert(cd_nelmts >= 4);
	unsigned char bytes[8];	
	*dimSize = cd_values[0];
	*dataType = cd_values[1];

	switch(*dimSize)
	{
	case 1:
		intToBytes_bigEndian(bytes, cd_values[2]);
		intToBytes_bigEndian(&bytes[4], cd_values[3]);
		if(sizeof(size_t)==4)
			*r1 = (unsigned int)bytesToLong_bigEndian(bytes);
		else
			*r1 = (unsigned long)bytesToLong_bigEndian(bytes);
		*r2 = *r3 = *r4 = *r5 = 0;
		break;
	case 2:
		*r3 = *r4 = *r5 = 0;
		*r2 = cd_values[3];
		*r1 = cd_values[2];
		break;
	case 3:
		*r4 = *r5 = 0;
		*r3 = cd_values[4];
		*r2 = cd_values[3];
		*r1 = cd_values[2];
		break;
	case 4:
		*r5 = 0;
		*r4 = cd_values[5];
		*r3 = cd_values[4];
		*r2 = cd_values[3];
		*r1 = cd_values[2];	
		break;
	default: 
		*r5 = cd_values[6];
		*r4 = cd_values[5];
		*r3 = cd_values[4];
		*r2 = cd_values[3];
		*r1 = cd_values[2];		
	}
	//printf("ROCCI_cdArrayToMetaData: r5=%zu, r4=%zu, r3=%zu, r2=%zu, r1=%zu\n", *r5, *r4, *r3, *r2, *r1);	
}

void H5Z_ROCCI_Init(char* cfgFile)
{
	ROCCI_Init(cfgFile);
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	if(r1==0)
	{
		dataLength = 0;
	}
	else if(r2==0)
	{
		dataLength = r1;
	}
	else if(r3==0)
	{
		dataLength = r1*r2;
	}
	else if(r4==0)
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0)
	{
		dataLength = r1*r2*r3*r4;
	}
	else
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

/**
 * @brief		check dimension and correct it if needed
 * @return 	0 (didn't change dimension)
 * 					1 (dimension is changed)
 * 					2 (dimension is problematic)
 **/
int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension)
{
	int dimensionCorrected = 0;
	int dim = computeDimension(r5, r4, r3, r2, r1);
	correctedDimension[0] = r1;
	correctedDimension[1] = r2;
	correctedDimension[2] = r3;
	correctedDimension[3] = r4;
	correctedDimension[4] = r5;
	size_t* c = correctedDimension;
	if(dim==1)
	{
		if(r1<1)
			return 2;
	}
	else if(dim==2)
	{
		if(r2==1)
		{
			c[1]= 0;
			dimensionCorrected = 1;
		}
		if(r1==1) //remove this dimension
		{
			c[0] = c[1];
			c[1] = c[2];
			dimensionCorrected = 1;
		}
	}
	else if(dim==3)
	{
		if(r3==1)
		{
			c[2] = 0;
			dimensionCorrected = 1;
		}
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			dimensionCorrected = 1;
		}
	}
	else if(dim==4)
	{
		if(r4==1)
		{
			c[3] = 0;
			dimensionCorrected = 1;
		}
		if(r3==1)
		{
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
	}
	else if(dim==5)
	{
		if(r5==1)
		{
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r4==1)
		{
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r3==1)
		{
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
	}

	return dimensionCorrected;

}

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
        int dimension;
        if(r1==0)
        {
                dimension = 0;
        }
        else if(r2==0)
        {
                dimension = 1;
        }
        else if(r3==0)
        {
                dimension = 2;
        }
        else if(r4==0)
        {
                dimension = 3;
        }
        else if(r5==0)
        {
                dimension = 4;
        }
        else
        {
                dimension = 5;
        }
        return dimension;
}


void ROCCI_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	unsigned char bytes[8] = {0};
	unsigned long size;	
	int dim = computeDimension(r5, r4, r3, r2, r1);
	cd_values[0] = dim;
	cd_values[1] = dataType;	//0: FLOAT ; 1: DOUBLE ; 2,3,4,....: INTEGER....
	switch(dim)
	{
	case 1:
		size = (unsigned long)r1;
		longToBytes_bigEndian(bytes, size);
		cd_values[2] = bytesToInt_bigEndian(bytes);
		cd_values[3] = bytesToInt_bigEndian(&bytes[4]);	
		*cd_nelmts = 4;
		break;
	case 2:
		cd_values[2] = (unsigned int) r2;
		cd_values[3] = (unsigned int) r1;	
		*cd_nelmts = 4;
		break;
	case 3:
		cd_values[2] = (unsigned int) r3;
		cd_values[3] = (unsigned int) r2;
		cd_values[4] = (unsigned int) r1;	
		*cd_nelmts = 5;
		break;
	case 4:
		cd_values[2] = (unsigned int) r4;	
		cd_values[3] = (unsigned int) r3;
		cd_values[4] = (unsigned int) r2;
		cd_values[5] = (unsigned int) r1;	
		*cd_nelmts = 6;
		break;
	default:
		cd_values[2] = (unsigned int) r5;		
		cd_values[3] = (unsigned int) r4;	
		cd_values[4] = (unsigned int) r3;
		cd_values[5] = (unsigned int) r2;
		cd_values[6] = (unsigned int) r1;
		*cd_nelmts = 7;	
	}	
}

/**
 * to be used in compression, and to be called outside H5Z_filter_rocci().
 * */
 
 void ROCCI_refreshDimForCdArray(int dataType, size_t old_cd_nelmts, unsigned int *old_cd_values, size_t* new_cd_nelmts, unsigned int **new_cd_values, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
 {
	 unsigned char bytes[8] = {0};
	 *new_cd_values = (unsigned int*)malloc(sizeof(unsigned int)*16);
	 memset(*new_cd_values, 0, sizeof(unsigned int)*16);
	 
	//correct dimension if needed
	size_t _r[5];
	filterDimension(r5, r4, r3, r2, r1, _r);
	size_t _r5 = _r[4];
	size_t _r4 = _r[3];
	size_t _r3 = _r[2];
	size_t _r2 = _r[1];
	size_t _r1 = _r[0];	
 
	int i =0;
	int oldDim = computeDimension(r5, r4, r3, r2, r1);
	int newDim = computeDimension(_r5, _r4, _r3, _r2, _r1);
	(*new_cd_values)[0] = newDim;
	(*new_cd_values)[1] = dataType;

	
	switch(newDim)
	{		
		case 1:
			longToBytes_bigEndian(bytes, (unsigned long)r1);
			(*new_cd_values)[2] = bytesToInt_bigEndian(bytes);
			(*new_cd_values)[3] = bytesToInt_bigEndian(&bytes[4]);	
			if(old_cd_nelmts==0)
				*new_cd_nelmts = 4;
			else
			{
				(*new_cd_values)[4] = old_cd_values[0];
				(*new_cd_values)[5] = old_cd_values[1];
				(*new_cd_values)[6] = old_cd_values[2];
				(*new_cd_values)[7] = old_cd_values[3];
				(*new_cd_values)[8] = old_cd_values[4];
				(*new_cd_values)[9] = old_cd_values[5];
				(*new_cd_values)[10] = old_cd_values[6];
				(*new_cd_values)[11] = old_cd_values[7];
				(*new_cd_values)[12] = old_cd_values[8];
				(*new_cd_values)[13] = old_cd_values[9];
				(*new_cd_values)[14] = old_cd_values[10];
				*new_cd_nelmts = 15;
			}
			break;		
		case 2:
			(*new_cd_values)[2] = (unsigned int) _r2;
			(*new_cd_values)[3] = (unsigned int) _r1;	
			if(old_cd_nelmts==0)
				*new_cd_nelmts = 4;
			else
			{
				(*new_cd_values)[4] = old_cd_values[0];
				(*new_cd_values)[5] = old_cd_values[1];
				(*new_cd_values)[6] = old_cd_values[2];
				(*new_cd_values)[7] = old_cd_values[3];
				(*new_cd_values)[8] = old_cd_values[4];
				(*new_cd_values)[9] = old_cd_values[5];
				(*new_cd_values)[10] = old_cd_values[6];
				(*new_cd_values)[11] = old_cd_values[7];
				(*new_cd_values)[12] = old_cd_values[8];
				(*new_cd_values)[13] = old_cd_values[9];
				(*new_cd_values)[14] = old_cd_values[10];
				*new_cd_nelmts = 15;
			}
			break;
		case 3:
			(*new_cd_values)[2] = (unsigned int) _r3;
			(*new_cd_values)[3] = (unsigned int) _r2;
			(*new_cd_values)[4] = (unsigned int) _r1;	
			if(old_cd_nelmts==0)
				*new_cd_nelmts = 5;
			else
			{
				(*new_cd_values)[5] = old_cd_values[0];
				(*new_cd_values)[6] = old_cd_values[1];
				(*new_cd_values)[7] = old_cd_values[2];
				(*new_cd_values)[8] = old_cd_values[3];
				(*new_cd_values)[9] = old_cd_values[4];
				(*new_cd_values)[10] = old_cd_values[5];
				(*new_cd_values)[11] = old_cd_values[6];
				(*new_cd_values)[12] = old_cd_values[7];
				(*new_cd_values)[13] = old_cd_values[8];
				(*new_cd_values)[14] = old_cd_values[9];
				(*new_cd_values)[15] = old_cd_values[10];				
				*new_cd_nelmts = 16;
			}
			break;
		case 4:
			(*new_cd_values)[2] = (unsigned int) _r4;	
			(*new_cd_values)[3] = (unsigned int) _r3;
			(*new_cd_values)[4] = (unsigned int) _r2;
			(*new_cd_values)[5] = (unsigned int) _r1;	
			if(old_cd_nelmts==0)
				*new_cd_nelmts = 6;
			else
			{
				(*new_cd_values)[6] = old_cd_values[0];
				(*new_cd_values)[7] = old_cd_values[1];
				(*new_cd_values)[8] = old_cd_values[2];
				(*new_cd_values)[9] = old_cd_values[3];
				(*new_cd_values)[10] = old_cd_values[4];
				(*new_cd_values)[11] = old_cd_values[5];
				(*new_cd_values)[12] = old_cd_values[6];
				(*new_cd_values)[13] = old_cd_values[7];
				(*new_cd_values)[14] = old_cd_values[8];
				(*new_cd_values)[15] = old_cd_values[9];
				(*new_cd_values)[16] = old_cd_values[10];				
				*new_cd_nelmts = 17;
				break;
			}
		default:
			(*new_cd_values)[2] = (unsigned int) _r5;		
			(*new_cd_values)[3] = (unsigned int) _r4;	
			(*new_cd_values)[4] = (unsigned int) _r3;
			(*new_cd_values)[5] = (unsigned int) _r2;
			(*new_cd_values)[6] = (unsigned int) _r1;
			if(old_cd_nelmts==0)
				*new_cd_nelmts = 7;
			else
			{
				(*new_cd_values)[7] = old_cd_values[0];
				(*new_cd_values)[8] = old_cd_values[1];
				(*new_cd_values)[9] = old_cd_values[2];
				(*new_cd_values)[10] = old_cd_values[3];
				(*new_cd_values)[11] = old_cd_values[4];						
				(*new_cd_values)[12] = old_cd_values[5];
				(*new_cd_values)[13] = old_cd_values[6];
				(*new_cd_values)[14] = old_cd_values[7];
				(*new_cd_values)[15] = old_cd_values[8];
				(*new_cd_values)[16] = old_cd_values[9];
				(*new_cd_values)[17] = old_cd_values[10];
				*new_cd_nelmts = 18;
			}
	}
 }

void ROCCI_errConfigToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int error_bound_mode, double abs_error, double rel_error, double pw_rel_error, double psnr, double fixedCR)
{
	*cd_values = (unsigned int*)malloc(sizeof(unsigned int)*11);
	int k = 0;
	(*cd_values)[k++] = error_bound_mode;
	unsigned char b[8];
	doubleToBytes(b, abs_error);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b+4);
	doubleToBytes(b, rel_error);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b);	
	(*cd_values)[k++] = bytesToInt32_bigEndian(b+4);
	doubleToBytes(b, pw_rel_error);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b);		
	(*cd_values)[k++] = bytesToInt32_bigEndian(b+4);
	doubleToBytes(b, psnr);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b+4);
	doubleToBytes(b, fixedCR);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b);
	(*cd_values)[k++] = bytesToInt32_bigEndian(b+4);
	*cd_nelmts = k;
}

static herr_t H5Z_rocci_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id)
{

	detectSysEndianType();
	setDataEndianType(LITTLE_ENDIAN_DATA);

	//printf("start in H5Z_rocci_set_local, dcpl_id = %d\n", dcpl_id);
	static char const *_funcname_ = "H5Z_rocci_set_local";
	size_t r5=0,r4=0,r3=0,r2=0,r1=0, dsize;
	int i, ndims, ndims_used = 0;	
	hsize_t dims[H5S_MAX_RANK], dims_used[5] = {0,0,0,0,0};	
	herr_t retval = 0;
	H5T_class_t dclass;
	H5T_sign_t dsign;
	unsigned int flags = 0;
	size_t mem_cd_nelmts = 11, cd_nelmts = 0;
	unsigned int mem_cd_values[22]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 

	//H5Z_FILTER_ROCCI
	//note that mem_cd_nelmts must be non-zero, otherwise, mem_cd_values cannot be filled.
	if (0 > H5Pget_filter_by_id(dcpl_id, H5Z_FILTER_ROCCI, &flags, &mem_cd_nelmts, mem_cd_values, 0, NULL, NULL))
		H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current ROCCI cd_values");	
	
	//the duplicated H5Pget_filter_by_id() is not a careless mistake....
	//The first H5Pget_filter_by_id() is to get the number of cd_values, and the second one is to fill the cd_values.
	//if (0 > H5Pget_filter_by_id(dcpl_id, H5Z_FILTER_ROCCI, &flags, &mem_cd_nelmts, mem_cd_values, 0, NULL, NULL))
	//	H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current ROCCI cd_values");

	//printf("22mem_cd_nelmts=%zu\n", mem_cd_nelmts);
	//for(int i=0;i<mem_cd_nelmts;i++)
	//	printf("22mem_cd_values[%d]: %d\n", i, mem_cd_values[i]);

	if(mem_cd_nelmts==0) //this means that the error information is missing from the cd_values
	{
		// printf("ERR INFO NOT FOUND IN CDVALS\n");
		H5Z_ROCCI_Init(cfgFile);
	}
	else //this means that the error information is included in the cd_values
	{
		ROCCI_Init(NULL);
		// printf("ERR INFO FOUND IN CDVALS\n");
		
	}

	herr_t ret = H5Zregister(H5Z_ROCCI); 
	if(ret < 0)
		printf("Error: H5Zregister(H5Z_ROCCI) faild.");
	
	int dataType = ROCCI_FLOAT;
	
	if (0 > (dclass = H5Tget_class(type_id)))
		H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");

	if (0 == (dsize = H5Tget_size(type_id)))
		H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");

	if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims, 0)))
		H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");
		
	for (i = 0; i < ndims; i++)
		dims_used[i] = dims[i];
	
	//printf("dclass=%d, H5T_FLOAT=%d, H5T_INTEGER=%d\n", dclass, H5T_FLOAT, H5T_INTEGER);
	if (dclass == H5T_FLOAT)
		dataType = dsize==4? ROCCI_FLOAT: ROCCI_DOUBLE;
	else if(dclass == H5T_INTEGER)
	{
		if (0 > (dsign = H5Tget_sign(type_id)))
			H5Z_ROCCI_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");		
		if(dsign == H5T_SGN_NONE) //unsigned
		{
			switch(dsize)
			{
			case 1:
				dataType = ROCCI_UINT8;
				break;
			case 2:
				dataType = ROCCI_UINT16;
				break;
			case 4:
				dataType = ROCCI_UINT32;
				break;
			case 8:
				dataType = ROCCI_UINT64;
				break;
			}
		}
		else
		{
			switch(dsize)
			{
			case 1:
				dataType = ROCCI_INT8;
				break;
			case 2:
				dataType = ROCCI_INT16;
				break;
			case 4:
				dataType = ROCCI_INT32;
				break;
			case 8:
				dataType = ROCCI_INT64;
				break;
			}			
		}
	}
	else
	{
		//printf("Error: dclass...\n");
		H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
	}
	
	unsigned int* cd_values = NULL;
	// printf("MEM %i\n", mem_cd_nelmts);
	if(mem_cd_nelmts!=0 && mem_cd_nelmts!=11)
	{
		H5Epush(H5E_DEFAULT,__FILE__, "H5Z_rocci_set_local", __LINE__, H5E_ERR_CLS, H5E_ARGS, H5E_BADVALUE, "Wrong number of cd_values: The new version has 9 integer elements in cd_values. Please check 'test/print_h5repack_args' to get the correct cd_values.");
		H5Eprint(H5E_DEFAULT, stderr);
		return -1;
	}
	ROCCI_refreshDimForCdArray(dataType, mem_cd_nelmts, mem_cd_values, &cd_nelmts, &cd_values, dims_used[4], dims_used[3], dims_used[2], dims_used[1], dims_used[0]);
	//printf("cd_nelmts=%zu\n", cd_nelmts);
	
//	ROCCI_metaDataToCdArray(&cd_nelmts, &cd_values, dataType, _r5, _r4, _r3, _r2, _r1);
	/* Now, update cd_values for the filter */
	if (0 > H5Pmodify_filter(dcpl_id, H5Z_FILTER_ROCCI, flags, cd_nelmts, cd_values))
		H5Z_ROCCI_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");	

	free(cd_values);

	retval = 1;
done:
	return retval;
}


int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[])
{
	int result = 0; //0 means no-error-information-in-cd_values; 1 means cd_values contains error information, -1 means wrong cd_values
	int dimSize = cd_values[0];
	//printf("nc_nelmts = %d\n", cd_nelmts);i
	//printf("dimSize=%d, cd_nelmts=%d, \n", dimSize, cd_nelmts);
	switch(dimSize)
	{
	case 1:
		if(cd_nelmts>4) //we are using two integers to store 1D data size length (long type)
			result = 1;
		break;
	case 2:
		if(cd_nelmts>4)
			result = 1;
		break;
	case 3:
		if(cd_nelmts>5)
			result = 1;
		break;
	case 4:
		if(cd_nelmts>6)
			result = 1;
		break;
	case 5:
		if(cd_nelmts>7)
			result = 1;
		break;
	}
	return result;
}

static size_t H5Z_filter_rocci(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf)
{
	//printf("start in H5Z_filter_rocci, cd_nelmts=%d\n", cd_nelmts);
	//H5Z_ROCCI_Init_Default();
	
	size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;
	int dimSize = 0, dataType = 0;

	if(cd_nelmts==0) //this is special data such as string, which should not be treated as values.
		return nbytes;
	
	int withErrInfo = checkCDValuesWithErrors(cd_nelmts, cd_values);

	int error_mode = 0;
	double abs_error = 0, rel_error = 0, pw_rel_error = 0, psnr = 0, ratio = 0;
	if(withErrInfo)
		ROCCI_cdArrayToMetaDataErr(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1, &error_mode, &abs_error, &rel_error, &pw_rel_error, &psnr, &ratio);
	else
		ROCCI_cdArrayToMetaData(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1);
	
	/*int i=0;
	for(i=0;i<cd_nelmts;i++)
		printf("cd_values[%d]=%u\n", i, cd_values[i]);
	printf("dimSize=%d, r1=%u, r2=%u, r3=%u, r4=%u, r5=%u\n", dimSize, r1, r2, r3, r4, r5);*/
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1); 

	if(nbEle < 20)
		return nbytes;

	if (flags & H5Z_FLAG_REVERSE) 
	{ 
		//cost_start();
		/* decompress data */
		if(dataType == ROCCI_FLOAT)//==0
		{
			float* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);

			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(float);
		}
		else if(dataType == ROCCI_DOUBLE)//==1
		{
			double* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(double);			
		}
		else if(dataType == ROCCI_INT8)
		{
			char* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(char);
		}
		else if(dataType == ROCCI_UINT8)
		{
			unsigned char* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned char);
		}
		else if(dataType == ROCCI_INT16)
		{
			short* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(short);
		}
		else if(dataType == ROCCI_UINT16)
		{
			unsigned short* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned short);
		}
		else if(dataType == ROCCI_INT32)
		{
			int* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(int);
		}
		else if(dataType == ROCCI_UINT32)
		{
			unsigned int* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned int);
		}
		else if(dataType == ROCCI_INT64)
		{
			long* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(long);
		}
		else if(dataType == ROCCI_UINT64)
		{
			unsigned long* data = ROCCI_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned long);
		}
		else
		{
			printf("Decompression error: unknown data type: %d\n", dataType);
			exit(0);
		}
		//cost_end();
		//printf("decompression time = %lf, decompression rate = %lf\n", totalCost, 1.0*nbEle*sizeof(float)/totalCost);
	}
	else //compression
	{
		size_t outSize = 0;
	
		//printf("r5=%d, r4=%d, r3=%d, r2=%d, r1=%d, dataType=%d\n", r5, r4, r3, r2, r1, dataType);
		//cost_start();
		if(dataType == ROCCI_FLOAT)//==0
		{
			float* data = (float*)(*buf);
			unsigned char *bytes = NULL;
			//printf("error=%f, withErrInfo=%d\n", confparams_cpr->relBoundRatio, withErrInfo);
			if(withErrInfo)
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_DOUBLE)//==1
		{
			double* data = (double*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				//printf("dataType=%d, error_mode=%d, abs_err=%f, rel_err=%f\n", dataType, error_mode, abs_error, rel_error);
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
			{
				//printf("....\n");
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			}
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;			
		}
		else if(dataType == ROCCI_INT8)
		{
			char* data = (char*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_UINT8)
		{
			unsigned char* data = (unsigned char*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_INT16)
		{
			short* data = (short*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_UINT16)
		{
			unsigned short* data = (unsigned short*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_INT32)
		{
			int* data = (int*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_UINT32)
		{
			unsigned int* data = (unsigned int*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_INT64)
		{
			long* data = (long*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == ROCCI_UINT64)
		{
			unsigned long* data = (unsigned long*)(*buf);
			unsigned char *bytes = NULL;
			if(withErrInfo)
			{
				bytes = ROCCI_compress_args(dataType, data, &outSize, error_mode, abs_error, rel_error, pw_rel_error, psnr, ratio, r5, r4, r3, r2, r1);
			}
			else
				bytes = ROCCI_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else 
		{
			printf("Compression error: unknown data type: %d\n", dataType);
			exit(0);
		}
		//cost_end();
		//printf("compression time = %lf, compression rate = %lf\n", totalCost, 1.0*nbEle*sizeof(float)/totalCost);
	}
	
	//H5Z_ROCCI_Finalize();
	return *buf_size;
}

void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	switch(dim)
	{
	case 1: 
		dims[0] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
			chunk[0] = r1;
		else
			chunk[0] = 2147483648;//2^31
		break;
	case 2:
		dims[0] = r2;
		dims[1] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r2;
			chunk[1] = r1;
		}
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}	
		break;
	case 3:
		dims[0] = r3;
		dims[1] = r2;
		dims[2] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r3;
			chunk[1] = r2;
			chunk[2] = r1;
		}		
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}
		break;
	case 4:
		dims[0] = r4;
		dims[1] = r3;
		dims[2] = r2;
		dims[3] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r4;
			chunk[1] = r3;
			chunk[2] = r2;
			chunk[3] = r1;
		}		
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}
		break;
	default:
		dims[0] = r5;
		dims[1] = r4;
		dims[2] = r3;
		dims[3] = r2;
		dims[4] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r5;
			chunk[1] = r4;
			chunk[2] = r3;
			chunk[3] = r2;
			chunk[4] = r1;
		}		
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}
	}
}




