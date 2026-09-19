#ifndef SOIL_CLI_H
#define SOIL_CLI_H

#include <limits.h>   /* for PATH_MAX fallback */
#include "soil_config.h"  /* brings NDISC, THETA_MIN, etc. */

/* -------- time-column mode for output files ------------------------ */
typedef enum {
    TC_INDEX = 0,   /* step counter */
    TC_JD,          /* strict Julian date */
    TC_YMDH         /* YYYY\u2011MM\u2011DD HH:MM */
} TimeColMode;

#ifndef MAX_PATH
#define MAX_PATH 1024
#endif

/* -------- driver / command line option bundle ---------------------- */
typedef struct {
    char forcing_path[MAX_PATH];
    char outdir[MAX_PATH];

    int  verbosity;          /* 0 silent, 1 MB/step, >1 verbose */
    int  write_theta;
    int  write_fluxes;
    int  write_mb;

    char solver[16];          /* dsbm or noahmp */
    int  apply_fc_perc_threshold; /* apply theta_fc threshold to bottom drainage */

    /* initialisation */
    int  use_init_wt_depth;
    double init_wt_depth_m;
    int  use_init_target_storage;
    double init_target_storage_m;

    /* lookup table */
    int  use_lut;
    int  lut_n;
    double lut_Theta_min;

    /* geometry overrides */
    int  have_dz;
    double dz_override[NDISC];
    int  have_zc;
    double zc_override[NDISC];

    /* output time column */
    TimeColMode timecol;

    /* flags set by parse_args() */
    int  have_forcing;
    int  have_outdir;
    int  have_cfe_config;
    int  tindex_start;                 /* --tindex-start */
    char  params_path[MAX_PATH];        /* optional separate params file */
    int   have_params;                  /* bool flag, like have_dz */
    
} DriverOpts;

/* parse_args implemented in bmi_soil_driver.c */
int parse_args(int argc, char **argv, DriverOpts *opt);

#endif /* SOIL_CLI_H */
