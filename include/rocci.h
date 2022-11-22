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

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif





#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ROCCI_H  ----- */
