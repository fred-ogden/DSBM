#define _POSIX_C_SOURCE 200809L

/* src/bmi_soil_driver.c
 *
 * Command-line driver for the stateless soil kernel.
 * Reads hourly forcing CSV records:
 *
 *   datetime,rainfall_mm_per_h,potential_et_mm_per_h
 *   YYYY/MM/DD hh:mm:ss,rainfall_mm_per_h,potential_et_mm_per_h
 *
 * Writes CSV time series for theta, fluxes, and volume balance,
 * plus a volume-balance summary file. Positive vertical flux is downward.
 *
 * Build example:
 *   cc -O3 -std=c11 -Wall -Wextra -Iinclude \
 *      src/soil_helpers.c src/soil_kernel_stateless.c src/bmi_soil_driver.c -o soil_driver
 *
 * Example run:
 *   ./soil_driver \
 *     --forcing forcing.txt \
 *     --outdir out \
 *     --verbosity 1 \
 *     --use-lut --lut-n 400 --lut-Theta-min 1e-6 \
 *     --init-wt-depth 1.0 \
 *     --write-theta --write-fluxes --write-volbal
 */

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
#include <time.h>
#include <math.h>

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
        "  --config <file>                Soil parameter configuration file.\n"
        "  --solver <dsbm|noahmp>          Soil solver (default dsbm).\n"
        "  --apply-fc-perc-threshold       Restrict bottom drainage to water above theta_fc.\n"
        "  --verbosity <int>              0 silent, 1 step volume balance, >1 detailed (default 1).\n"
        "  --write-theta                  Write theta_timeseries.csv (comma delimited).\n"
        "  --write-fluxes                 Write fluxes_timeseries.csv (comma delimited).\n"
        "  --write-volbal                 Write volbal_timeseries.csv and volbal_summary.out.\n"
        "\n"
        "Timestamp format:\n"
        "  --timestamp <timestep|datetime|juliandate>\n"
        "                                 Output timestamp format (default datetime).\n"
        "      timestep   : integer timestep index\n"
        "      datetime   : YYYY-MM-DD hh:mm:ss\n"
        "      juliandate : astronomical Julian Date (JD; noon-based origin)\n"
        "  --timestep-start <int>          Starting index for --timestamp timestep (default 0).\n"
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
        "Forcing file format (CSV):\n"
        "  datetime,rainfall_mm_per_h,potential_et_mm_per_h\n"
        "  YYYY/MM/DD hh:mm:ss,rainfall_mm_per_h,potential_et_mm_per_h\n"
        "  The header row is required. Lines beginning with '#' are ignored.\n"
        "\n"
        "Notes:\n"
        "  * Outputs are comma-delimited CSV with a header line starting with '#'.\n"
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
    o->timecol = TC_DATETIME;   // NEW
    strcpy(o->solver, "dsbm");
    o->timestep_start = 0;           // NEW
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
        case TC_DATETIME:
            fprintf(fp, "%04d-%02d-%02d %02d:%02d:00", YYYY, MM, DD, HH, Mi);
            break;
        case TC_JULIANDATE: {
            double jd = julian_date_ymdhm(YYYY, MM, DD, HH, Mi);
            fprintf(fp, "%.10f", jd);
            break;
        }
        case TC_TIMESTEP:
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
        } else if (!strcmp(a, "--config") && i+1 < argc) {
            strncpy(o->config_path, argv[++i], MAX_PATH-1);
            o->have_config = 1;
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
        } else if (!strcmp(a, "--write-volbal")) {
            o->write_volbal = 1;
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
        } else if (strcmp(a, "--timestamp") == 0 && i+1 < argc) {
            const char *m = argv[++i];
            if      (!strcasecmp(m, "timestep"))   o->timecol = TC_TIMESTEP;
            else if (!strcasecmp(m, "datetime"))   o->timecol = TC_DATETIME;
            else if (!strcasecmp(m, "juliandate")) o->timecol = TC_JULIANDATE;
            else {
                fprintf(stderr, "Unknown --timestamp value '%s'\n", m);
                return -1;
            }
        } else if (!strcmp(a, "--timestep-start") && i+1 < argc) {
            if (parse_int(argv[++i], &o->timestep_start)) return -1;
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
    if (!o->use_init_wt_depth && !o->use_init_target_storage && !o->have_config) {
        fprintf(stderr, "Initialization is required on the command line or in --config.\n");
        return -1;
    }
    return 0;
}

// ---------------------- Parameters & geometry ----------------------

typedef struct {
    int have_soil_depth_m;
    double soil_depth_m;
    int calculate_phi_sat_from_ksat;
    int have_field_capacity_pressure_ratio;
    double field_capacity_pressure_ratio;
    int have_initial_storage_m;
    double initial_storage_m;
} SoilConfigOptions;

static int load_config_file(const char *path, SoilParameters *par, SoilConfigOptions *config)
{
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    char line[512];
    long line_number = 0;
    memset(config, 0, sizeof(*config));

    while (fgets(line, sizeof(line), fp)) {
        line_number++;
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '#') continue;

        char key[160], value[160];
        if (sscanf(p, " %159[^= \t] = %159[^ \t\n]", key, value) != 2) {
            fprintf(stderr, "ERROR: malformed soil config record at line %ld:\n%s",
                    line_number, line);
            fclose(fp);
            return -1;
        }

        char *units = strchr(value, '[');
        if (units) *units = '\0';

        if (!strcmp(key, "soil_saturated_capillary_head_calc_from_ksat")) {
            if (!strcmp(value, "TRUE")) config->calculate_phi_sat_from_ksat = 1;
            else if (!strcmp(value, "FALSE")) config->calculate_phi_sat_from_ksat = 0;
            else {
                fprintf(stderr, "ERROR: %s must be TRUE or FALSE at line %ld.\n",
                        key, line_number);
                fclose(fp);
                return -1;
            }
            continue;
        }

        char *endptr = NULL;
        double val = strtod(value, &endptr);
        if (endptr == value || *endptr != '\0') {
            fprintf(stderr, "ERROR: invalid value for %s at line %ld: %s\n",
                    key, line_number, value);
            fclose(fp);
            return -1;
        }

        if (!strcmp(key, "soil_depth_m")) {
            config->soil_depth_m = val; config->have_soil_depth_m = 1;
        } else if (!strcmp(key, "soil_Clapp_Hornberger_exponent_b")) {
            par->b_exp = val;
        } else if (!strcmp(key, "soil_sat_hydraulic_conductivity_cm_per_h")) {
            par->K_sat_cm_per_h = val;
        } else if (!strcmp(key, "soil_sat_capillary_head_cm")) {
            par->phi_sat_cm = val;
        } else if (!strcmp(key, "soil_effective_porosity")) {
            par->theta_sat = val;
        } else if (!strcmp(key, "soil_field_capacity_Pcap_over_Patm_0_1")) {
            config->field_capacity_pressure_ratio = val;
            config->have_field_capacity_pressure_ratio = 1;
        } else if (!strcmp(key, "soil_reservoir_rate_const_to_subsurface_lateral_flow")) {
            par->klf_m_per_h = val;
        } else if (!strcmp(key, "soil_to_gw_percolation_rate_limiter_0_to_1")) {
            par->perc_limiter_0_to_1 = val;
        } else if (!strcmp(key, "state_soil_reservoir_init_storage_m")) {
            config->initial_storage_m = val; config->have_initial_storage_m = 1;
        } else {
            fprintf(stderr, "ERROR: unknown soil config parameter at line %ld: %s\n",
                    line_number, key);
            fclose(fp);
            return -1;
        }
    }
    fclose(fp);

    if (config->calculate_phi_sat_from_ksat)
        par->phi_sat_cm = 10.415 * pow(par->K_sat_cm_per_h, -0.3266);

    if (config->have_field_capacity_pressure_ratio) {
        const double atmospheric_pressure_head_cm = 1033.2274528;
        double field_capacity_head_cm =
            config->field_capacity_pressure_ratio * atmospheric_pressure_head_cm;
        par->theta_fc = (field_capacity_head_cm <= par->phi_sat_cm)
                      ? par->theta_sat
                      : par->theta_sat *
                        pow(par->phi_sat_cm / field_capacity_head_cm, 1.0 / par->b_exp);
        par->theta_aet_eq_pet = par->theta_fc;
    }
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
    if (opt->timecol == TC_DATETIME) fprintf(fp, "# datetime");
    else if (opt->timecol == TC_JULIANDATE) fprintf(fp, "# juliandate");
    else fprintf(fp, "# timestep");

    for (int i = 0; i < NDISC; i++) fprintf(fp, ",theta_%d", i+1);
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
    if (opt->timecol == TC_DATETIME) fprintf(fp, "# datetime");
    else if (opt->timecol == TC_JULIANDATE) fprintf(fp, "# juliandate");
    else fprintf(fp, "# timestep");

    fprintf(fp, ",rain_mm_per_h,rain_into_soil_m,percolation_m");
    for (int i = 0; i < NDISC; i++) fprintf(fp, ",AET_disc_%d_m", i+1);
    for (int i = 0; i < NDISC; i++) fprintf(fp, ",lateral_disc_%d_m", i+1);
    for (int i = 0; i < NDISC; i++) fprintf(fp, ",interface_rate_%d_m_per_h", i+1);
    for (int i = 0; i < NDISC; i++) fprintf(fp, ",interface_vol_%d_m", i+1);
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

static void write_volbal_header(FILE *fp, const DriverOpts *opt)
{
    if (opt->timecol == TC_DATETIME) fprintf(fp, "# datetime");
    else if (opt->timecol == TC_JULIANDATE) fprintf(fp, "# juliandate");
    else fprintf(fp, "# timestep");
    fprintf(fp, ",in_rain_m,excess_m,perc_m,AET_m,lateral_m,delta_storage_m,residual_m\n");
}

// COMMENTED OUT .csv OUTPUT
//static void write_volbal_header(FILE *fp)
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

    SoilConfigOptions config;
    memset(&config, 0, sizeof(config));
    if (opt.have_config) {
        if (load_config_file(opt.config_path, &par, &config))
            die("failed to read soil config file");
    }

    fill_geometry_from_opts(&opt, &geom);
    double soil_depth_m = soil_depth_from_dz(&geom);
    if (config.have_soil_depth_m && fabs(config.soil_depth_m - soil_depth_m) > 1.0e-10) {
        fprintf(stderr, "ERROR: soil_depth_m=%.10g does not match sum(dz)=%.10g m.\n",
                config.soil_depth_m, soil_depth_m);
        return EXIT_FAILURE;
    }

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
    } else if (opt.use_init_target_storage || config.have_initial_storage_m) {
        double target_storage_m = opt.use_init_target_storage
                                ? opt.init_target_storage_m
                                : config.initial_storage_m;
        double zwt;
        initialize_hydrostatic_from_storage(soil_depth_m,
                                            par.theta_sat, par.phi_sat_cm, par.b_exp,
                                            geom.zc, target_storage_m, &zwt, theta);
        if (opt.verbosity > 0)
            fprintf(stdout, "Initialized hydrostatic z_wt = %.6f m from target storage %.10g m\n",
                    zwt, target_storage_m);
    } else {
        die("no soil initialization specified");
    }

    // Compute initial psi/K
    compute_props_with_option_stateless(opt.use_lut ? &lut : NULL,
                                        theta, psi_init, K_init,
                                        par.theta_r, par.theta_sat,
                                        par.K_sat_cm_per_h, par.phi_sat_cm, par.b_exp,
                                        hint);

    double initial_soil_volume_m = 0.0;
    for (int i = 0; i < NDISC; i++) {
        initial_soil_volume_m += theta[i] * geom.dz[i];
    }

    // Outputs
    char path_theta[MAX_PATH], path_flux[MAX_PATH], path_volbal[MAX_PATH], path_sum[MAX_PATH];
    FILE *fp_theta = NULL, *fp_flux = NULL, *fp_volbal = NULL, *fp_sum = NULL;

    if (opt.write_theta) {
        fp_theta = open_csv(opt.outdir, "theta_timeseries.csv", path_theta);
        write_theta_header(fp_theta, &opt);
    }
    if (opt.write_fluxes) {
        fp_flux = open_csv(opt.outdir, "fluxes_timeseries.csv", path_flux);
        write_fluxes_header(fp_flux, &opt);
    }
    if (opt.write_volbal) {
        fp_volbal = open_csv(opt.outdir, "volbal_timeseries.csv", path_volbal);
        write_volbal_header(fp_volbal, &opt);
        fp_sum = open_csv(opt.outdir, "volbal_summary.out", path_sum);
    }


    // Forcing file
    FILE *fp_forcing = fopen(opt.forcing_path, "r");
    if (!fp_forcing) die_errno(opt.forcing_path);

    // Accumulators for simulation volume balance
    double cum_input_rain_m = 0.0;
    double cum_infiltrated_rain_m = 0.0;
    double cum_excess_m  = 0.0;
    double cum_perc_m    = 0.0;
    double cum_AET_m     = 0.0;
    double cum_lat_m     = 0.0;

    // State structs reused per step
    SoilStateIn  sin = {0};
    SoilStateOut sout = {0};
    SoilForcing  forcing = {0};
    SoilFluxes   flux = {0};
    TimestepSoilVolumeBalance volbal = {0};

    // this is used when writing output using a timestep index
    // EXTRA? long step_idx = opt.timestep_start;

    // fill constant parts
    // sin.theta_in will be filled each step from theta[]
    // sin.psi_in, sin.K_in filled from psi_init/K_init before the step
    // sin.ch_lut_hint_in from hint[]
    // copy initial hints
    for (int i = 0; i < NDISC; i++) {
        sout.ch_lut_hint_out[i] = hint[i];
    }

    // Main loop over forcing lines

/* ---------- Main loop over forcing CSV records ---------- */

    char line[1024];
    long forcing_line_number = 0;
    long n_steps = 0;
    long step_idx = opt.timestep_start;
    long long total_subtimesteps = 0;
    struct timespec wall_start, wall_end;
    if (clock_gettime(CLOCK_MONOTONIC, &wall_start) != 0) die_errno("clock_gettime");

    if (fgets(line, sizeof(line), fp_forcing) == NULL) {
        fprintf(stderr, "ERROR: forcing file is empty.\n");
        fclose(fp_forcing);
        return 1;
    }
    forcing_line_number++;


    const char *expected_header =
        "datetime,rainfall_mm_per_h,potential_et_mm_per_h";

    char *pita = line;
    while (*pita == ' ' || *pita == '\t') pita++;

    if (strncmp(pita, expected_header, strlen(expected_header)) != 0) {
        fprintf(stderr,
                "ERROR: forcing file line %ld: expected CSV header:\n"
                "%s\n",
                forcing_line_number, expected_header);
        fclose(fp_forcing);
        return 1;
    }


    while (fgets(line, sizeof(line), fp_forcing)) {
        forcing_line_number++;

        // Skip comments and blank lines.
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '#') continue;

        int YYYY, MM, DD, HH, Mi, SS;
        double P_mm, PET_mm;
        char extra;
        int nread = sscanf(p, "%d/%d/%d %d:%d:%d,%lf,%lf %c",
                           &YYYY, &MM, &DD, &HH, &Mi, &SS,
                           &P_mm, &PET_mm, &extra);
        if (nread != 8) {
            fprintf(stderr,
                    "ERROR: malformed forcing record at line %ld:\n%s",
                    forcing_line_number, line);
            fclose(fp_forcing);
            return 1;
        }
        if (SS < 0 || SS > 59) {
            fprintf(stderr, "ERROR: invalid seconds at forcing line %ld: %d\n",
                    forcing_line_number, SS);
            fclose(fp_forcing);
            return 1;
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
                &sin, &forcing, &sout, &flux, &volbal, NULL);
            if (rc != 0) die("noahmp_soil_step_one_hour_stateless failed");
        } else {
            rc = soil_step_one_hour_stateless(
                &ctrl, &geom, &par, (opt.use_lut ? &lut : NULL),
                &sin, &forcing, &sout, &flux, &volbal, NULL);
            if (rc != 0) die("soil_step_one_hour_stateless failed");
        }

        total_subtimesteps += flux.n_sub_used;

        // Update state for next step
        for (int i = 0; i < NDISC; i++) theta[i] = sout.theta_out[i];
        for (int i = 0; i < NDISC; i++) hint[i]  = sout.ch_lut_hint_out[i];

        compute_props_with_option_stateless(opt.use_lut ? &lut : NULL,
            theta, psi_init, K_init,
            par.theta_r, par.theta_sat,
            par.K_sat_cm_per_h, par.phi_sat_cm, par.b_exp,
            hint);

        // Accumulate simulation volume-balance terms
        cum_input_rain_m += forcing.rain_mm_per_h / 1000.0;
        cum_infiltrated_rain_m += volbal.in_rain_m;
        cum_excess_m  += volbal.excess_m;
        cum_perc_m    += volbal.perc_m;
        cum_AET_m     += volbal.AET_m;
        cum_lat_m     += volbal.lateral_m;

        // Verbose stdout
        if (opt.verbosity >= 1) {
            fprintf(stdout,
                "%04d-%02d-%02d %02d:%02d  rain_in=%.6e  perc=%.6e  AET=%.6e  lat=%.6e  dS=%.6e  res=%.3e  nsub=%d\n",
                YYYY, MM, DD, HH, Mi,
                volbal.in_rain_m, volbal.perc_m, volbal.AET_m, volbal.lateral_m, volbal.delta_storage_m,
                volbal.residual_m, flux.n_sub_used);
        }
        if (opt.verbosity > 1) {
            double th_min = theta[0], th_max = theta[0];
            for (int i = 1; i < NDISC; i++) { if (theta[i] < th_min) th_min = theta[i]; if (theta[i] > th_max) th_max = theta[i]; }
            fprintf(stdout, "   theta_min=%.6f theta_max=%.6f\n", th_min, th_max);
        }

        // ---- Outputs (space-delimited) ----
        if (fp_theta) {
            fprint_time_col(fp_theta, &opt, step_idx, YYYY, MM, DD, HH, Mi);
            for (int i = 0; i < NDISC; i++) fprintf(fp_theta, ",%.10g", theta[i]);
            fprintf(fp_theta, "\n");
        }

        if (fp_flux) {
            fprint_time_col(fp_flux, &opt, step_idx, YYYY, MM, DD, HH, Mi);
            fprintf(fp_flux, ",%.10g,%.10g,%.10g",
                    forcing.rain_mm_per_h, flux.rain_into_soil_m, flux.percolation_to_gw_m);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, ",%.10g", flux.AET_by_disc_m[i]);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, ",%.10g", flux.lateral_by_disc_m[i]);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, ",%.10g", flux.interface_rate_m_per_h[i]);
            for (int i = 0; i < NDISC; i++) fprintf(fp_flux, ",%.10g", flux.interface_vol_m[i]);
            fprintf(fp_flux, "\n");
        }

        if (fp_volbal) {
            fprint_time_col(fp_volbal, &opt, step_idx, YYYY, MM, DD, HH, Mi);
            fprintf(fp_volbal, " %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
                    volbal.in_rain_m, volbal.excess_m, volbal.perc_m, volbal.AET_m, volbal.lateral_m,
                    volbal.delta_storage_m, volbal.residual_m);
        }

        step_idx++;
        n_steps++;
    } // <-- make sure this closes the while loop


    if (clock_gettime(CLOCK_MONOTONIC, &wall_end) != 0) die_errno("clock_gettime");
    double wall_clock_seconds =
        (double)(wall_end.tv_sec - wall_start.tv_sec) +
        1.0e-9 * (double)(wall_end.tv_nsec - wall_start.tv_nsec);

    fclose(fp_forcing);

    double final_soil_volume_m = 0.0;
    for (int i = 0; i < NDISC; i++) {
        final_soil_volume_m += theta[i] * geom.dz[i];
    }

    double volume_balance_residual_m =
        initial_soil_volume_m + cum_input_rain_m
        - cum_excess_m - cum_AET_m - cum_lat_m - cum_perc_m
        - final_soil_volume_m;

    if (fp_sum) {
        const char *solver_name = !strcmp(opt.solver, "noahmp") ? "Noah-MP" : "DSBM";
        double mean_subtimesteps = (n_steps > 0)
            ? (double)total_subtimesteps / (double)n_steps : 0.0;

        fprintf(fp_sum, "Solver: %s\n\n", solver_name);
        fprintf(fp_sum, "**************** %s SOIL VOLUME BALANCE ****************\n\n", solver_name);
        fprintf(fp_sum, " Initial soil volume                    = %10.4f m\n", initial_soil_volume_m);
        fprintf(fp_sum, " Precipitation input (QINSUR)            = %10.4f m\n", cum_input_rain_m);
        fprintf(fp_sum, " Surface precipitation excess            = %10.4f m\n", cum_excess_m);
        fprintf(fp_sum, " Infiltrated precipitation               = %10.4f m\n", cum_infiltrated_rain_m);
        fprintf(fp_sum, " Evapotranspiration                      = %10.4f m\n", cum_AET_m);
        fprintf(fp_sum, " Lateral flow generated                  = %10.4f m\n", cum_lat_m);
        fprintf(fp_sum, " Percolation to g.w.                     = %10.4f m\n", cum_perc_m);
        fprintf(fp_sum, " Final soil volume                      = %10.4f m\n", final_soil_volume_m);
        fprintf(fp_sum, " Volume balance residual                = %10.4e m\n", volume_balance_residual_m);
        if (cum_input_rain_m != 0.0) {
            fprintf(fp_sum, " Residual as percent of input rainfall  = %10.4e percent\n",
                    100.0 * volume_balance_residual_m / cum_input_rain_m);
        } else {
            fprintf(fp_sum, " Residual as percent of input rainfall  = NaN percent (no input rainfall)\n");
        }

        fprintf(fp_sum, "\n***************** SIMULATION INFORMATION ******************\n\n");
        fprintf(fp_sum, " Simulation duration                    = %10ld h\n", n_steps);
        fprintf(fp_sum, " Output timestep                        = hourly\n");
        fprintf(fp_sum, " Total number of sub-timesteps          = %10lld\n", total_subtimesteps);
        fprintf(fp_sum, " Mean sub-timesteps per output step     = %10.3f\n", mean_subtimesteps);
        fprintf(fp_sum, " Simulation wall-clock time             = %10.3f s\n", wall_clock_seconds);
        fclose(fp_sum);
    }
    if (fp_theta) fclose(fp_theta);
    if (fp_flux)  fclose(fp_flux);
    if (fp_volbal) fclose(fp_volbal);

    if (opt.use_lut) soil_free_ch_lut(&lut);

    return 0;
}
