/*
 * util.h
 *
 * Public interface for general driver utility functions.
 *
 * Discrete Soil Moisture Balance Model (DSBM)
 *
 * Author:
 *   Fred L. Ogden, Ph.D., P.E.
 *   NOAA/National Weather Service
 *
 * This software was developed by an employee of the United States
 * Government as part of official duties and is not subject to
 * copyright protection in the United States under 17 U.S.C. Section 105.
 *
 * License:
 *   Apache License, Version 2.0
 *   SPDX-License-Identifier: Apache-2.0
 *
 * See the repository LICENSE file for additional information.
 */

#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>
#include "soil_cli.h"   /* needs TimeColMode */

int  util_ensure_dir(const char *path);

double util_to_julian(int Y,int M,int D,int HH,int Mi,int Sec);

void fprintf_print_time_col(FILE       *fp,
                            TimeColMode mode,
                            long        idx,
                            int         Y,int M,int D,
                            int         HH,int Mi);

#endif
