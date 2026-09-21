/* src/util.c */
#include "util.h"
#include <sys/stat.h>
#include <errno.h>

#ifdef _WIN32
#include <direct.h>   /* _mkdir */
#endif

/* mkdir -p equivalent */
int util_ensure_dir(const char *path)
{
    struct stat st;
    if (stat(path,&st)==0) return S_ISDIR(st.st_mode) ? 0 : -1;
#ifdef _WIN32
    return _mkdir(path);
#else
    return mkdir(path,0777);
#endif
}

/* Gregorian to JD (days) */
double util_to_julian(int Y,int M,int D,int HH,int Mi,int Sec)
{
    int A=(14-M)/12, Yp=Y+4800-A, Mp=M+12*A-3;
    long J=D+(153*Mp+2)/5+365L*Yp+Yp/4-Yp/100+Yp/400-32045;
    return J + (HH-12)/24.0 + Mi/1440.0 + Sec/86400.0;
}

/* time-column printer */
void fprintf_time_col(FILE *fp, TimeColMode mode, long idx,
                      int Y,int M,int D,int HH,int Mi)
{
    switch(mode) {
    case TC_TIMESTEP: fprintf(fp,"%ld", idx); break;
    case TC_JULIANDATE:    fprintf(fp,"%.8f", util_to_julian(Y,M,D,HH,Mi,0)); break;
    case TC_DATETIME:  fprintf(fp,"%04d-%02d-%02d %02d:%02d", Y,M,D,HH,Mi); break;
    }
}

