/**
 *  @file H5Z_ROCCI.h
 *  @author Sheng Di
 *  @date July, 2022
 *  @brief Header file for the H5Z_ROCCI.c.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _H5Z_ROCCI
#define _H5Z_ROCCI

#include <stdio.h>
#include <hdf5.h>
#include "rocci.h"

#define H5Z_FILTER_ROCCI 32017
#define MAX_CHUNK_SIZE 4294967295 //2^32-1
static hid_t H5Z_ROCCI_ERRCLASS = -1;

#ifdef __cplusplus
extern "C" {
#endif


extern int load_conffile_flag;
extern int init_rocci_flag;

extern char cfgFile[256];

/* convenience macro to handle errors */
#define ERROR(FNAME)                                              \
do {                                                              \
    int _errno = errno;                                           \
    fprintf(stderr, #FNAME " failed at line %d, errno=%d (%s)\n", \
        __LINE__, _errno, _errno?strerror(_errno):"ok");          \
    return 1;                                                     \
} while(0)

#define H5Z_ROCCI_PUSH_AND_GOTO(MAJ, MIN, RET, MSG)     \
do                                                    \
{                                                     \
	H5Epush(H5E_DEFAULT,__FILE__,_funcname_,__LINE__, \
		H5Z_ROCCI_ERRCLASS,MAJ,MIN,MSG);                \
	retval = RET;                                     \
	goto done;                                        \
} while(0)


void ROCCI_refreshDimForCdArray(int dataType, size_t old_cd_nelmts, unsigned int *old_cd_values, size_t* new_cd_nelmts, unsigned int **new_cd_values, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void ROCCI_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1);
void ROCCI_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void ROCCI_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
int* error_bound_mode, double* abs_error, double* rel_error, double* pw_rel_error, double* psnr);

void ROCCI_errConfigToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int error_bound_mode, double abs_error, double rel_error, double pw_rel_error, double psnr, double fixedCR);

int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[]);
static size_t H5Z_filter_rocci(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf);
static herr_t H5Z_rocci_set_local(hid_t dcpl_id, hid_t type_id, hid_t space_id);


void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _H5Z_ROCCI_metadata  ----- */
