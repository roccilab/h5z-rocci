/**
 *  @file rocci_ByteToolkit.h
 *  @author Sheng Di
 *  @date July, 2023
 *  @brief Header file for the rocci_ByteToolkit.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _rocci_ByteToolkit_H
#define _rocci_ByteToolkit_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "rocci_defines.h"

//ByteToolkit.c

int sysEndianType;
int dataEndianType;

typedef union lint16 {
    unsigned short usvalue;
    short svalue;
    unsigned char byte[2];
} lint16;

typedef union lint32 {
    int ivalue;
    unsigned int uivalue;
    unsigned char byte[4];
} lint32;

typedef union lint64 {
    int64_t lvalue;
    uint64_t ulvalue;
    unsigned char byte[8];
} lint64;

typedef union ldouble {
    double value;
    uint64_t lvalue;
    unsigned char byte[8];
} ldouble;

typedef union lfloat {
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
    uint16_t int16[2];
} lfloat;

void detectSysEndianType();

void symTransform_4bytes(unsigned char data[4]);
void symTransform_8bytes(unsigned char data[8]);

void setSysEndianType(int endianType);
void setDataEndianType(int endianType);

extern unsigned short bytesToUInt16_bigEndian(unsigned char* bytes);
extern unsigned int bytesToUInt32_bigEndian(unsigned char* bytes);
extern uint64_t bytesToUInt64_bigEndian(unsigned char* b);

extern short bytesToInt16_bigEndian(unsigned char* bytes);
extern int bytesToInt32_bigEndian(unsigned char* bytes);
extern int64_t bytesToInt64_bigEndian(unsigned char* b);
extern int bytesToInt_bigEndian(unsigned char* bytes);

extern void intToBytes_bigEndian(unsigned char *b, unsigned int num);
int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension);
extern void int64ToBytes_bigEndian(unsigned char *b, uint64_t num);
extern void int32ToBytes_bigEndian(unsigned char *b, uint32_t num);
extern void int16ToBytes_bigEndian(unsigned char *b, uint16_t num);

extern int64_t bytesToLong_bigEndian(unsigned char* b);
extern void longToBytes_bigEndian(unsigned char *b, uint64_t num);
int64_t doubleToOSEndianLong(double value);
int floatToOSEndianInt(float value);
extern short getExponent_float(float value);
extern short getPrecisionReqLength_float(float precision);
extern short getExponent_double(double value);
extern short getPrecisionReqLength_double(double precision);
unsigned char numberOfLeadingZeros_Int(int i);
unsigned char numberOfLeadingZeros_Long(int64_t i);
unsigned char getLeadingNumbers_Int(int v1, int v2);
unsigned char getLeadingNumbers_Long(int64_t v1, int64_t v2);
short bytesToShort(unsigned char* bytes);
void shortToBytes(unsigned char* b, short value);
int bytesToInt(unsigned char* bytes);
int64_t bytesToLong(unsigned char* bytes);
extern float bytesToFloat(unsigned char* bytes);
extern void floatToBytes(unsigned char *b, float num);
extern double bytesToDouble(unsigned char* bytes);
extern void doubleToBytes(unsigned char *b, double num);
int getMaskRightCode(int m);
extern int getLeftMovingCode(int kMod8);
extern int getRightMovingSteps(int kMod8, int resiBitLength);
extern int getRightMovingCode(int kMod8, int resiBitLength);
short* convertByteDataToShortArray(unsigned char* bytes, size_t byteLength);
unsigned short* convertByteDataToUShortArray(unsigned char* bytes, size_t byteLength);

void convertShortArrayToBytes(short* states, size_t stateLength, unsigned char* bytes);
void convertUShortArrayToBytes(unsigned short* states, size_t stateLength, unsigned char* bytes);
void convertIntArrayToBytes(int* states, size_t stateLength, unsigned char* bytes);
void convertUIntArrayToBytes(unsigned int* states, size_t stateLength, unsigned char* bytes);
void convertLongArrayToBytes(int64_t* states, size_t stateLength, unsigned char* bytes);
void convertULongArrayToBytes(uint64_t* states, size_t stateLength, unsigned char* bytes);

extern size_t bytesToSize(unsigned char* bytes);
extern void sizeToBytes(unsigned char* outBytes, size_t size);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ByteToolkit_H  ----- */

