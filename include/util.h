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
