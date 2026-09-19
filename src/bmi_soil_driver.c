// src/bmi_soil_driver.c
//
// Command-line driver for the stateless soil kernel.
// Reads hourly forcing lines:
//
//   YYYY MM DD HH Mi  Precip_mm  PET_mm
//
// Writes CSV time series for theta, fluxes, and mass balance,
// plus a summary text file. Positive vertical flux is downward.
//
// Build example:
//   cc -O3 -std=c11 -Wall -Wextra -Iinclude \
//      src/soil_helpers.c src/soil_kernel_stateless.c src/bmi_soil_driver.c -o soil_driver
//
// Example run:
//   ./soil_driver \
//     --forcing forcing.txt \
//     --outdir out \
//     --verbosity 1 \
//     --use-lut --lut-n 400 --lut-Theta-min 1e-6 \
//     --init-wt-depth 1.0 \
//     --write-theta --write-fluxes --write-mb
//
// Required headers from this project:
#include <strings.h>
#include "soil_data_types.h"
#include "soil_helpers.h"
#include "soil_kernel_stateless.h"
#include "noahmp_soilwater_stateless.h"
#include "soil_cli.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>

#ifndef MAX_PATH
#define MAX_PATH 4096
#endif


// --------------------------- Utilities -----------------------------

static void die(const char *msg)
{
    fprintf(stderr, "FATAL: %s\n", msg);
    exit(EXIT_FAILURE);
}

static void die_errno(const char *msg)
{
    fprintf(stderr, "FATAL: %s : %s\n", msg, strerror(errno));
    exit(EXIT_FAILURE);
}

static int ensure_outdir(const char *path)
{
#if defined(_WIN32)
    int rc = _mkdir(path);
    if (rc == 0 || errno == EEXIST) return 0;
    return -1;
#else
    int rc = mkdir(path, 0775);
    if (rc == 0 || errno == EEXIST) return 0;
    return -1;
#endif
}


static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s --forcing <file> --outdir <dir> [options]\n"
        "\n"
        "Options:\n"
        "  --params <file>                Parameters text file (key value per line).\n"
        "  --solver <dsbm|noahmp>          Soil solver (default dsbm).\n"
        "  --apply-fc-perc-threshold       Restrict bottom drainage to water above theta_fc.\n"
        "  --verbosity <int>              0 silent, 1 step MB, >1 detailed (default 1).\n"
        "  --write-theta                  Write theta_timeseries.csv (space delimited).\n"
        "  --write-fluxes                 Write fluxes_timeseries.csv (space delimited).\n"
        "  --write-mb                     Write massbal_timeseries.csv and summary.\n"
        "\n"
        "Time column (choose format):\n"
        "  --time-col <calendar|jd|index> Time column mode (default calendar).\n"
        "      calendar : prints YYYY MM DD HH Mi\n"
        "      jd       : prints strict Julian Date (JD; noon-based origin)\n"
        "      index    : prints integer step index\n"
        "  --tindex-start <int>           Starting index when --time-col index (default 0).\n"
        "  --utc-offset-hours <h>         Hours to subtract from local clock to get UTC\n"
        "                                 before JD conversion (default 0). Use if input\n"
        "                                 timestamps are not UTC.\n"
        "\n"
        "Initialization (choose one):\n"
        "  --init-wt-depth <m>            Hydrostatic profile given water table depth.\n"
        "  --init-target-storage <m>      Hydrostatic profile matching target storage.\n"
        "\n"
        "Lookup table (optional):\n"
        "  --use-lut                      Use CH LUT (faster than analytic).\n"
        "  --lut-n <int>                  Number of LUT points (default 400).\n"
        "  --lut-Theta-min <val>          Minimum Theta for LUT grid (default 1e-6).\n"
        "\n"
        "Geometry (optional if compiled in):\n"
        "  --dz v1,v2,...,vN              Disc thicknesses (m) for NDISC discs.\n"
        "  --zc v1,v2,...,vN              Disc center depths (m).\n"
        "\n"
        "Forcing file format (space or general whitespace; comments allowed):\n"
        "  Lines beginning with '#' are ignored.\n"
        "  YYYY MM DD HH Mi  Precip_mm  PET_mm   (Mi typically 0)\n"
        "\n"
        "Notes:\n"
        "  * Outputs are space-delimited with a header line starting with '#'.\n"
        "  * JD uses astronomical convention (day changes at noon). If your input\n"
        "    times are local, set --utc-offset-hours accordingly to produce strict JD.\n",
        prog);
}


static int parse_csv_list(const char *s, double *out, int n)
{
    // Parse comma-separated list of n doubles.
    int count = 0;
    const char *p = s;
    while (*p && count < n) {
        char *end = NULL;
        double v = strtod(p, &end);
        if (end == p) return -1;
        out[count++] = v;

        p = end;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == ',') p++;
        while (*p == ' ' || *p == '\t') p++;
    }
    return (count == n) ? 0 : -1;
}

static int parse_int(const char *s, int *v_out)
{
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (end == s || *end != '\0') return -1;
    *v_out = (int)v;
    return 0;
}

static int parse_double(const char *s, double *v_out)
{
    char *end = NULL;
    double v = strtod(s, &end);
    if (end == s || *end != '\0') return -1;
    *v_out = v;
    return 0;
}

static void default_opts(DriverOpts *o)
{
    memset(o, 0, sizeof(*o));
    o->verbosity = 1;
    o->lut_n = 400;
    o->lut_Theta_min = 1e-6;
    o->timecol = TC_YMDH;   // NEW
    strcpy(o->solver, "dsbm");
    o->tindex_start = 0;           // NEW
    o->forcing_path[0] = '\0';
    o->outdir[0] = '\0';
}

static double julian_date_ymdhm(int Y, int M, int D, int H, int Mi)
/* Strict astronomical JD. Assumes input time is UTC.
   JD day changes at noon: add (H-12)/24 + Mi/1440. */
{
    int a = (14 - M)/12;
    int y = Y + 4800 - a;
    int m = M + 12*a - 3;
    long JDN = D + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045; // at 0h UT
    double frac = (H - 12)/24.0 + Mi/1440.0;
    return (double)JDN + frac;
}

static void fprint_time_col(FILE *fp, const DriverOpts *opt,
                            long step_index, int YYYY, int MM, int DD, int HH, int Mi)
{
    switch (opt->timecol) {
        case TC_YMDH:
            fprintf(fp, "%04d %02d %02d %02d %02d", YYYY, MM, DD, HH, Mi);
            break;
        case TC_JD: {
            double jd = julian_date_ymdhm(YYYY, MM, DD, HH, Mi);
            fprintf(fp, "%.10f", jd);
            break;
        }
        case TC_INDEX:
            fprintf(fp, "%ld", step_index);
            break;
    }
}

int parse_args(int argc, char **argv, DriverOpts *o)
{
    default_opts(o);
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--forcing") && i+1 < argc) {
            strncpy(o->forcing_path, argv[++i], MAX_PATH-1);
        } else if (!strcmp(a, "--outdir") && i+1 < argc) {
            strncpy(o->outdir, argv[++i], MAX_PATH-1);
        } else if (!strcmp(a, "--params") && i+1 < argc) {
            strncpy(o->params_path, argv[++i], MAX_PATH-1);
            o->have_params = 1;
        } else if (!strcmp(a, "--solver") && i+1 < argc) {
            const char *solver = argv[++i];
            if (strcmp(solver, "dsbm") && strcmp(solver, "noahmp")) return -1;
            strncpy(o->solver, solver, sizeof(o->solver)-1);
            o->solver[sizeof(o->solver)-1] = '\0';
        } else if (!strcmp(a, "--apply-fc-perc-threshold")) {
            o->apply_fc_perc_threshold = 1;
        } else if (!strcmp(a, "--verbosity") && i+1 < argc) {
            if (parse_int(argv[++i], &o->verbosity)) return -1;
        } else if (!strcmp(a, "--write-theta")) {
            o->write_theta = 1;
        } else if (!strcmp(a, "--write-fluxes")) {
            o->write_fluxes = 1;
        } else if (!strcmp(a, "--write-mb")) {
            o->write_mb = 1;
        } else if (!strcmp(a, "--init-wt-depth") && i+1 < argc) {
            if (parse_double(argv[++i], &o->init_wt_depth_m)) return -1;
            o->use_init_wt_depth = 1;
        } else if (!strcmp(a, "--init-target-storage") && i+1 < argc) {
            if (parse_double(argv[++i], &o->init_target_storage_m)) return -1;
            o->use_init_target_storage = 1;
        } else if (!strcmp(a, "--use-lut")) {
            o->use_lut = 1;
        } else if (!strcmp(a, "--lut-n") && i+1 < argc) {
            if (parse_int(argv[++i], &o->lut_n)) return -1;
        } else if (!strcmp(a, "--lut-Theta-min") && i+1 < argc) {
            if (parse_double(argv[++i], &o->lut_Theta_min)) return -1;
        } else if (!strcmp(a, "--dz") && i+1 < argc) {
            if (parse_csv_list(argv[++i], o->dz_override, NDISC)) return -1;
            o->have_dz = 1;
        } else if (!strcmp(a, "--zc") && i+1 < argc) {
            if (parse_csv_list(argv[++i], o->zc_override, NDISC)) return -1;
            o->have_zc = 1;
        } else if (strcmp(a,"--timecol")==0 && i+1 < argc) {
            const char *m = argv[++i];
            if      (!strcasecmp(m,"index")) o->timecol = TC_INDEX;
            else if (!strcasecmp(m,"jd"))    o->timecol = TC_JD;
            else if (!strcasecmp(m,"ymdh"))  o->timecol = TC_YMDH;
            else {
                fprintf(stderr,"Unknown --timecol value '%s'\n", m);
                return -1;
            }
        } else if (!strcmp(a, "--tindex-start") && i+1 < argc) {
            if (parse_int(argv[++i], &o->tindex_start)) return -1;
        } else {
            fprintf(stderr, "Unknown or incomplete option: %s\n", a);
            return -1;
        }
    }

    if (o->forcing_path[0] == '\0' || o->outdir[0] == '\0') return -1;
    if (o->use_init_wt_depth && o->use_init_target_storage) {
        fprintf(stderr, "Choose only one of --init-wt-depth or --init-target-storage.\n");
        return -1;
    }
    if (!o->use_init_wt_depth && !o->use_init_target_storage) {
        fprintf(stderr, "One initialization option is required.\n");
        return -1;
    }
    return 0;
}

// ---------------------- Parameters & geometry ----------------------

static int load_params_file(const char *path, SoilParameters *par)
{
    // Simple key value format; lines like:
    // theta_r 0.02
    // theta_sat 0.45
    // theta_fc 0.30
    // theta_wp 0.10
    // theta_aet_eq_pet 0.30
    // K_sat_cm_per_h 10
    // phi_sat_cm 10
    // b_exp 4.9
    // perc_limiter_0_to_1 1.0
    // klf_m_per_h 0.0
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    char key[128];
    double val;
    while (fscanf(fp, "%127s %lf", key, &val) == 2) {
        if (!strcmp(key, "theta_r")) par->theta_r = val;
        else if (!strcmp(key, "theta_sat")) par->theta_sat = val;
        else if (!strcmp(key, "theta_fc")) par->theta_fc = val;
        else if (!strcmp(key, "theta_wp")) par->theta_wp = val;
        else if (!strcmp(key, "theta_aet_eq_pet")) par->theta_aet_eq_pet = val;
        else if (!strcmp(key, "K_sat_cm_per_h")) par->K_sat_cm_per_h = val;
        else if (!strcmp(key, "phi_sat_cm")) par->phi_sat_cm = val;
        else if (!strcmp(key, "b_exp")) par->b_exp = val;
        else if (!strcmp(key, "perc_limiter_0_to_1")) par->perc_limiter_0_to_1 = val;
        else if (!strcmp(key, "klf_m_per_h")) par->klf_m_per_h = val;
        // ignore unknown keys
    }
    fclose(fp);
    return 0;
}

static void default_params(SoilParameters *par)
{
    memset(par, 0, sizeof(*par));
    par->theta_r = 0.0;
    par->theta_sat = 0.45;
    par->theta_fc = 0.30;
    par->theta_wp = 0.10;
    par->theta_aet_eq_pet = 0.30;
    par->K_sat_cm_per_h = 10.0;
    par->phi_sat_cm = 10.0;
    par->b_exp = 4.9;
    par->perc_limiter_0_to_1 = 1.0;
    par->klf_m_per_h = 0.0;
}

static void default_control(SoilControl *ctrl)
{
    memset(ctrl, 0, sizeof(*ctrl));
    ctrl->ndisc = NDISC;
    ctrl->deepest_root_disc = NDISC; // all discs in root zone by default
    ctrl->use_ch_lookup_table = 0;   // analytic by default
    ctrl->apply_fc_perc_threshold = 0; // native free drainage unless requested
    ctrl->dt_hours = 1.0;
}

static void fill_geometry_from_opts(const DriverOpts *o, SoilGeometry *geom)
{
    // If provided via CLI, copy; else zero-init (user must compile in constants).
    static double dz_local[NDISC];
    static double zc_local[NDISC];

    if (o->have_dz) {
        for (int i = 0; i < NDISC; i++) dz_local[i] = o->dz_override[i];
    } else {
        // Provide a safe default (uniform 0.25 m for 4 discs as an example)
        // Replace with your compiled geometry as needed.
        for (int i = 0; i < NDISC; i++) dz_local[i] = 1.0 / (double)NDISC;
    }

    if (o->have_zc) {
        for (int i = 0; i < NDISC; i++) zc_local[i] = o->zc_override[i];
    } else {
        double z = 0.0;
        for (int i = 0; i < NDISC; i++) {
            z += 0.5 * dz_local[i];
            zc_local[i] = z;
            z += 0.5 * dz_local[i];
        }
    }

    // Assign to const fields by casting away const here (safe for our process lifetime).
    memcpy((double*)geom->dz, dz_local, sizeof(dz_local));
    memcpy((double*)geom->zc, zc_local, sizeof(zc_local));
}

static double soil_depth_from_dz(const SoilGeometry *geom)
{
    double d = 0.0;
    for (int i = 0; i < NDISC; i++) d += geom->dz[i];
    return d;
}

// ---------------------------- I/O helpers --------------------------

static FILE* open_csv(const char *outdir, const char *name, char *fullpath)
{
    snprintf(fullpath, MAX_PATH, "%s/%s", outdir, name);
    FILE *fp = fopen(fullpath, "w");
    if (!fp) die_errno(fullpath);
    return fp;
}

static void write_theta_header(FILE *fp, const DriverOpts *opt)
{
    if (opt->timecol == TC_YMDH) fprintf(fp, "# YYYY MM DD HH Mi");
    else if (opt->timecol == TC_JD) fprintf(fp, "# JD");
    else fprintf(fp, "# step_index");

    for (int i = 0; i < NDISC; i++) fprintf(fp, " theta_%d", i+1);
    fprintf(fp, "\n");
}

// COMMENTED OUT .csv OUTPUT
//static void write_theta_header(FILE *fp)
//{
//    fprintf(fp, "YYYY,MM,DD,HH,Mi");
//    for (int i = 0; i < NDISC; i++) fprintf(fp, ",theta_%d", i+1);
//    fprintf(fp, "\n");
//}

static void write_fluxes_header(FILE *fp, const DriverOpts *opt)
{
    if (opt->timecol == TC_YMDH) fprintf(fp, "# YYYY MM DD HH Mi");
    else if (opt->timecol == TC_JD) fprintf(fp, "# JD");
    else fprintf(fp, "# step_index");

    fprintf(fp, " rain_mm_per_h rain_into_soil_m percolation_m");
    for (int i = 0; i < NDISC; i++) fprintf(fp, " AET_disc_%d_m", i+1);
    for (int i = 0; i < NDISC; i++) fprintf(fp, " lateral_disc_%d_m", i+1);
    for (int i = 0; i < NDISC; i++) fprintf(fp, " interface_rate_%d_m_per_h", i+1);
    for (int i = 0; i < NDISC; i++) fprintf(fp, " interface_vol_%d_m", i+1);
    fprintf(fp, "\n");
}

// COMMENTED OUT .csv OUTPUT
//static void write_fluxes_header(FILE *fp)
//{
//    fprintf(fp, "YYYY,MM,DD,HH,Mi");
//    fprintf(fp, ",rain_mm_per_h,rain_into_soil_m,percolation_m");
//    for (int i = 0; i < NDISC; i++) fprintf(fp, ",AET_disc_%d_m", i+1);
//    for (int i = 0; i < NDISC; i++) fprintf(fp, ",lateral_disc_%d_m", i+1);
//    for (int i = 0; i < NDISC; i++) fprintf(fp, ",interface_rate_%d_m_per_h", i+1);
//    for (int i = 0; i < NDISC; i++) fprintf(fp, ",interface_vol_%d_m", i+1);
//    fprintf(fp, "\n");
//}

static void write_mb_header(FILE *fp, const DriverOpts *opt)
{
    if (opt->timecol == TC_YMDH) fprintf(fp, "# YYYY MM DD HH Mi");
    else if (opt->timecol == TC_JD) fprintf(fp, "# JD");
    else fprintf(fp, "# step_index");
    fprintf(fp, " in_rain_m excess_m perc_m AET_m lateral_m delta_storage_m residual_m\n");
}

// COMMENTED OUT .csv OUTPUT
//static void write_mb_header(FILE *fp)
//{
//    fprintf(fp, "#YYYY,MM,DD,HH,Mi,in_rain_m,excess_m,perc_m,AET_m,lateral_m,delta_storage_m,residual_m\n");
//}

// ---------------------------- Main driver --------------------------

int main(int argc, char **argv)
{
    DriverOpts opt;
    if (parse_args(argc, argv, &opt)) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (ensure_outdir(opt.outdir)) die_errno("creating outdir");

    // Set up control, params, geometry
    SoilControl ctrl;   default_control(&ctrl);
    SoilParameters par; default_params(&par);
    SoilGeometry geom;  memset((void*)&geom, 0, sizeof(geom));

    ctrl.use_ch_lookup_table = opt.use_lut ? 1 : 0;
    ctrl.apply_fc_perc_threshold = opt.apply_fc_perc_threshold ? 1 : 0;

    if (opt.have_params) {
        if (load_params_file(opt.params_path, &par)) die("failed to read params file");
    }

    fill_geometry_from_opts(&opt, &geom);
    double soil_depth_m = soil_depth_from_dz(&geom);

    // Optional LUT
    SoilLookupTables lut = {0};
    if (opt.use_lut) {
        if (soil_build_ch_lut(opt.lut_n, opt.lut_Theta_min,
                              par.theta_r, par.theta_sat, par.b_exp,
                              par.phi_sat_cm, par.K_sat_cm_per_h, &lut) != 0)
        {
            die("Failed to build LUT");
        }
    }

    // Initialize theta
    double theta[NDISC];
    double psi_init[NDISC], K_init[NDISC];
    int    hint[NDISC]; for (int i = 0; i < NDISC; i++) hint[i] = -1;

    if (opt.use_init_wt_depth) {
        // Hydrostatic from specified water table depth
        double z_wt = opt.init_wt_depth_m;
        for (int i = 0; i < NDISC; i++) {
            double z = geom.zc[i];
            double phi_sat_m = par.phi_sat_cm / 100.0;
            double z_cf_top = z_wt - phi_sat_m;
            if (z >= z_cf_top) theta[i] = par.theta_sat;
            else {
                double psi_m = z_wt - z; if (psi_m < 0.0) psi_m = 0.0;
                theta[i] = theta_from_psi(psi_m, par.theta_sat, par.phi_sat_cm, par.b_exp);
            }
        }
    } else if (opt.use_init_target_storage) {
        double zwt;
        initialize_hydrostatic_from_storage(soil_depth_m,
                                            par.theta_sat, par.phi_sat_cm, par.b_exp,
                                            geom.zc,
                                            opt.init_target_storage_m,
                                            &zwt, theta);
        if (opt.verbosity > 0) {
            fprintf(stdout, "Initialized hydrostatic z_wt = %.6f m\n", zwt);
        }
    }

    // Compute initial psi/K
    compute_props_with_option_stateless(opt.use_lut ? &lut : NULL,
                                        theta, psi_init, K_init,
                                        par.theta_r, par.theta_sat,
                                        par.K_sat_cm_per_h, par.phi_sat_cm, par.b_exp,
                                        hint);

    // Outputs
    char path_theta[MAX_PATH], path_flux[MAX_PATH], path_mb[MAX_PATH], path_sum[MAX_PATH];
    FILE *fp_theta = NULL, *fp_flux = NULL, *fp_mb = NULL, *fp_sum = NULL;

    if (opt.write_theta) {
        fp_theta = open_csv(opt.outdir, "theta_timeseries.csv", path_theta);
        write_theta_header(fp_theta, &opt);
    }
    if (opt.write_fluxes) {
        fp_flux = open_csv(opt.outdir, "fluxes_timeseries.csv", path_flux);
        write_fluxes_header(fp_flux, &opt);
    }
    if (opt.write_mb) {
        fp_mb = open_csv(opt.outdir, "massbal_timeseries.csv", path_mb);
        write_mb_header(fp_mb, &opt);
        fp_sum = open_csv(opt.outdir, "massbal_summary.txt", path_sum);
    }


    // Forcing file
    FILE *fp_forcing = fopen(opt.forcing_path, "r");
    if (!fp_forcing) die_errno(opt.forcing_path);

    // Accumulators for global mass balance
    double cum_in_rain_m = 0.0;   // infiltrated
    double cum_excess_m  = 0.0;
    double cum_perc_m    = 0.0;
    double cum_AET_m     = 0.0;
    double cum_lat_m     = 0.0;
    double cum_resid_m   = 0.0;

    // State structs reused per step
    SoilStateIn  sin = {0};
    SoilStateOut sout = {0};
    SoilForcing  forcing = {0};
    SoilFluxes   flux = {0};
    TimestepSoilMassbal mb = {0};

    // this is used when writing output using a timestep index
    // EXTRA? long step_idx = opt.tindex_start;

    // fill constant parts
    // sin.theta_in will be filled each step from theta[]
    // sin.psi_in, sin.K_in filled from psi_init/K_init before the step
    // sin.ch_lut_hint_in from hint[]
    // copy initial hints
    for (int i = 0; i < NDISC; i++) {
        sout.ch_lut_hint_out[i] = hint[i];
    }

    // Main loop over forcing lines

/* ---------- Main loop over forcing lines (space delimited, '#' comments) ---------- */

    char line[1024];
    long n_steps = 0;
    long step_idx = opt.tindex_start;

    while (fgets(line, sizeof(line), fp_forcing)) {
        // Skip comments/blank
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '#') continue;

        int YYYY, MM, DD, HH, Mi;
        double P_mm, PET_mm;
        int nread = sscanf(p, "%d %d %d %d %d %lf %lf",
                           &YYYY, &MM, &DD, &HH, &Mi, &P_mm, &PET_mm);
        if (nread != 7) {
            if (opt.verbosity > 0) fprintf(stderr, "Warning: skip malformed line: %s", line);
            continue;
        }

        forcing.rain_mm_per_h = P_mm;
        forcing.pet_mm_per_h  = PET_mm;

        // Prepare input state
        for (int i = 0; i < NDISC; i++) {
            sin.theta_in[i]      = theta[i];
            sin.psi_in[i]        = psi_init[i];
            sin.K_in[i]          = K_init[i];
            sin.ch_lut_hint_in[i]= sout.ch_lut_hint_out[i];
        }

        // Run selected soil kernel
        int rc;
        if (!strcmp(opt.solver, "noahmp")) {
            rc = noahmp_soil_step_one_hour_stateless(
                &ctrl, &geom, &par,
                &sin, &forcing, &sout, &flux, &mb, NULL);
            if (rc != 0) die("noahmp_soil_step_one_hour_stateless failed");
        } else {
            rc = soil_step_one_hour_stateless(
                &ctrl, &geom, &par, (opt.use_lut ? &lut : NULL),
                &sin, &forcing, &sout, &flux, &mb, NULL);
            if (rc != 0) die("soil_step_one_hour_stateless failed");
        }

        // Update state for next step
        for (int i = 0; i < NDISC; i++) theta[i] = sout.theta_out[i];
        for (int i = 0; i < NDISC; i++) hint[i]  = sout.ch_lut_hint_out[i];

        compute_props_with_option_stateless(opt.use_lut ? &lut : NULL,
            theta, psi_init, K_init,
            par.theta_r, par.theta_sat,
            par.K_sat_cm_per_h, par.phi_sat_cm, par.b_exp,
            hint);

        // Global mass balances
        cum_in_rain_m += mb.in_rain_m;
        cum_excess_m  += mb.excess_m;
        cum_perc_m    += mb.perc_m;
        cum_AET_m     += mb.AET_m;
        cum_lat_m     += mb.lateral_m;
        cum_resid_m   += mb.residual_m;

        // Verbose stdout
        if (opt.verbosity >= 1) {
            fprintf(stdout,
                "%04d-%02d-%02d %02d:%02d  rain_in=%.6e  perc=%.6e  AET=%.6e  lat=%.6e  dS=%.6e  res=%.3e  nsub=%d\n",
                YYYY, MM, DD, HH, Mi,
                mb.in_rain_m, mb.perc_m, mb.AET_m, mb.lateral_m, mb.delta_storage_m,
                mb.residual_m, flux.n_sub_used);
        }
        if (opt.verbosity > 1) {
            double th_min = theta[0], th_max = theta[0];
            for (int i = 1; i < NDISC; i++) { if (theta[i] < th_min) th_min = theta[i]; if (theta[i] > th_max) th_max = theta[i]; }
            fprintf(stdout, "   theta_min=%.6f theta_max=%.6f\n", th_min, th_max);
        }

        // ---- Outputs (space-delimited) ----
        if (fp_theta) {
            fprint_time_col(fp_theta, &opt, step_idx, YYYY, MM, DD, HH, Mi);
            for (int i = 0; i < NDISC; i++) fprintf(fp_theta, " %.10g", theta[i]);
            fprintf(fp_theta, "\n");
        }

        if (fp_flux) {
            fprint_time_col(fp_flux, &opt, step_idx, YYYY, MM, DD, HH, Mi);
            fprintf(fp_flux, " %.10g %.10g %.10g",
                    forcing.rain_mm_per_h, flux.rain_into_soil_m, flux.percolation_to_gw_m);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, " %.10g", flux.AET_by_disc_m[i]);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, " %.10g", flux.lateral_by_disc_m[i]);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, " %.10g", flux.interface_rate_m_per_h[i]);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, " %.10g", flux.interface_vol_m[i]);
            fprintf(fp_flux, "\n");
        }

        if (fp_mb) {
            fprint_time_col(fp_mb, &opt, step_idx, YYYY, MM, DD, HH, Mi);
            fprintf(fp_mb, " %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
                    mb.in_rain_m, mb.excess_m, mb.perc_m, mb.AET_m, mb.lateral_m,
                    mb.delta_storage_m, mb.residual_m);
        }

        step_idx++;
        n_steps++;
    } // <-- make sure this closes the while loop


    fclose(fp_forcing);

    if (fp_sum) {
        fprintf(fp_sum, "Steps: %ld\n", n_steps);
        fprintf(fp_sum, "Cumulative infiltrated rain (m): %.10g\n", cum_in_rain_m);
        fprintf(fp_sum, "Cumulative excess (m): %.10g\n", cum_excess_m);
        fprintf(fp_sum, "Cumulative percolation (m): %.10g\n", cum_perc_m);
        fprintf(fp_sum, "Cumulative AET (m): %.10g\n", cum_AET_m);
        fprintf(fp_sum, "Cumulative lateral (m): %.10g\n", cum_lat_m);
        fprintf(fp_sum, "Cumulative residual (m): %.10g\n", cum_resid_m);
        fclose(fp_sum);
    }
    if (fp_theta) fclose(fp_theta);
    if (fp_flux)  fclose(fp_flux);
    if (fp_mb)    fclose(fp_mb);

    if (opt.use_lut) soil_free_ch_lut(&lut);

    return 0;
}
