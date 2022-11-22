/**
 *  @file rocci_defines.h
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief Header file for the dataCompression.c.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ROCCI_DEFINES_H
#define _ROCCI_DEFINES_H

#define ROCCI_VERNUM 0x0200
#define ROCCI_VER_MAJOR 1
#define ROCCI_VER_MINOR 0
#define ROCCI_VER_BUILD 0
#define ROCCI_VER_REVISION 0

#define ABS 0
#define REL 1
#define VR_REL 1  //alternative name to REL
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PSNR 4
#define NORM 5
#define FIX_RATE 6

#define PW_REL 10
#define ABS_AND_PW_REL 11
#define ABS_OR_PW_REL 12
#define REL_AND_PW_REL 13
#define REL_OR_PW_REL 14


#define ROCCI_FLOAT 0
#define ROCCI_DOUBLE 1
#define ROCCI_UINT8 2
#define ROCCI_INT8 3
#define ROCCI_UINT16 4
#define ROCCI_INT16 5
#define ROCCI_UINT32 6
#define ROCCI_INT32 7
#define ROCCI_UINT64 8
#define ROCCI_INT64 9

#define LITTLE_ENDIAN_DATA 0 //refers to the endian type of the data read from the disk
#define BIG_ENDIAN_DATA 1 //big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)

#define LITTLE_ENDIAN_SYSTEM 0 //refers to the endian type of the system
#define BIG_ENDIAN_SYSTEM 1


#define ROCCI_NO_BLOCK_FAST_CMPR 1
#define ROCCI_WITH_BLOCK_FAST_CMPR 2
#define ROCCI_RANDOMACCESS_FAST_CMPR 3
#define ROCCI_OPENMP_FAST_CMPR 4

//SUCCESS returning status
#define ROCCI_SCES 0  //successful
#define ROCCI_NSCS -1 //Not successful
#define ROCCI_FERR -2 //Failed to open input file
#define ROCCI_TERR -3 //wrong data type (should be only float or double)
#define ROCCI_DERR -4 //dimension error
#define ROCCI_MERR -5 //sz_mode error
#define ROCCI_BERR -6 //bound-mode error (should be only ABS, REL, ABS_AND_REL, ABS_OR_REL, or PW_REL)

#endif /* _ROCCI_DEFINES_H */
