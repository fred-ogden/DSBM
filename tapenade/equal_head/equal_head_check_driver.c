/*
 * equal_head_check_driver.c
 *
 * E11 (FLO 11/06/2026): is Tapenade's derivative of the equal-head limit
 * on downward transfers correct at every limiter activation that occurs
 * on the forcing record?
 *
 * The limiter solves for the downward transfer V at which the two discs
 * reach equal total head,
 *
 *   h(V, p) = ln(psi_upper(V)) - ln(psi_lower(V) + dz_int) = 0,
 *   psi_upper = phi_sat (theta_upper_after / theta_sat)^(-b),
 *   psi_lower = phi_sat (theta_lower_after / theta_sat)^(-b),
 *   theta_upper_after = theta_upper - V / dz_upper,
 *   theta_lower_after = theta_lower + V / dz_lower,
 *   dz_int = (dz_upper + dz_lower) / 2,
 *
 * with a safeguarded Newton iteration of fixed length (30 iterations).
 * By the implicit function theorem the exact derivative of the solution
 * with respect to any input p is
 *
 *   dV/dp = -(dh/dp) / (dh/dV),
 *
 * and dV/d(transfer_requested) = 0 (the requested transfer only brackets
 * the root).  Tapenade differentiates the iterations themselves, which
 * gives this derivative only if the iteration ends with Newton steps; the
 * census showed that bisection steps are common in the last iterations.
 *
 * Method:
 *   1. Run DSBM over the forcing record (analytic Clapp-Hornberger, fixed
 *      4 substeps, exponential lateral flow, no field-capacity
 *      percolation threshold) for one parameter set (k_lf, perc_limiter,
 *      Ksat, b), with the other soil properties as in
 *      configs/soil_params.dat and phi_sat, theta_fc and theta_wp derived
 *      as in the driver and in E10.  The kernel object is compiled with
 *      limit_downward_transfer_to_equal_head renamed to the recording
 *      wrapper below, so every limiter call with its arguments is
 *      recorded when the limit reduces the transfer.  The production
 *      kernel source is not changed.
 *   2. For every recorded activation that took the Newton path
 *      (0 < limited transfer < requested transfer), compare for each of
 *      the six active inputs:
 *        AD        Tapenade tangent (one call per input, unit seed)
 *        analytic  -(dh/dp)/(dh/dV) at the limited transfer
 *        FD        central finite difference of the production limiter
 *   3. For activations that set the transfer to zero (the discs were
 *      already at or past equal head), the exact derivative is zero for
 *      every input; check that the AD derivative is zero.
 *
 * Activations at which theta_upper_after or theta_lower_after sits at
 * a clamp of psi_from_theta() (the dry-end limit
 * CH_SCALED_SOIL_MOISTURE_MIN theta_sat, or theta_sat) are counted
 * separately: the analytic formula above assumes neither clamp is active.
 *
 * Usage:
 *   equal_head_check  output_directory  forcing_csv_file  k_lf_m_per_h  perc_limiter  K_sat_cm_per_h  b_exp
 *
 * Writes output_directory/activations.csv (one line per Newton-path
 * activation: inputs, head residual, and AD, analytic and FD derivatives)
 * and prints a summary.  Exit status 0 if every Newton-path activation
 * has AD within AD_ACCEPT_RELATIVE_ERROR of analytic for every input
 * and every zero-transfer activation has zero AD derivative; 2 otherwise.
 *
 * Terminology: soil moisture; discs.  ASCII only.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "soil_data_types.h"
#include "soil_helpers.h"
#include "dsbm_soilmoisture_stateless.h"
#include "equal_head_check.h"

/* ---- model setup, as configs/soil_params.dat and E10 ---- */
#define SOIL_THETA_SAT                     0.439
#define SOIL_DEPTH_M                       2.0
#define INITIAL_STORAGE_M                  0.7843498367003
#define FIELD_CAPACITY_PRESSURE_RATIO      0.333
#define WILTING_POINT_PRESSURE_RATIO       15.0
#define ATMOSPHERIC_PRESSURE_HEAD_CM       1033.2274528
#define PHI_SAT_FROM_KSAT_COEFFICIENT_CM   10.415
#define PHI_SAT_FROM_KSAT_EXPONENT         (-0.3266)

#define MAX_HOURS                          80000
#define MAX_ACTIVATIONS                    2000000

/* acceptance and finite-difference settings */
#define AD_ACCEPT_RELATIVE_ERROR           1.0e-8   /* AD versus analytic */
#define AD_ACCEPT_ABSOLUTE_DERIVATIVE      1.0e-12  /* floor for |analytic| in relative error */
#define FD_RELATIVE_STEP                   1.0e-6


/* ---- recorded limiter activations ---- */
typedef struct {
    long   hour_index;                          /* 0-based hour of the record */
    double limiter_inputs[N_LIMITER_INPUTS];    /* arguments, see equal_head_check.h */
    double limited_transfer_m;                  /* returned transfer (m) */
} LimiterActivation;

static LimiterActivation *activations = NULL;
static long n_activations = 0;
static long n_activations_dropped = 0;
static long n_limiter_calls = 0;
static long current_hour_index = 0;


/*
 * Recording wrapper.  The DSBM kernel object used by this program is
 * compiled with -Dlimit_downward_transfer_to_equal_head=equal_head_recording_limit,
 * so the kernel calls this function, which calls the production limiter
 * (from libsoil.a) and records the call when the limit reduced the
 * transfer.  The returned value is the production limiter's value.
 */
double equal_head_recording_limit(double transfer_requested_m,
                                  double theta_upper, double theta_lower,
                                  double dz_upper_m, double dz_lower_m,
                                  double theta_sat, double phi_sat_cm, double b_exp)
{
    double limited_transfer_m;
    LimiterActivation *record;

    limited_transfer_m = limit_downward_transfer_to_equal_head(
        transfer_requested_m, theta_upper, theta_lower,
        dz_upper_m, dz_lower_m, theta_sat, phi_sat_cm, b_exp);
    n_limiter_calls++;

    if (limited_transfer_m < transfer_requested_m) {
        if (n_activations < MAX_ACTIVATIONS) {
            record = &activations[n_activations];
            record->hour_index = current_hour_index;
            record->limiter_inputs[INPUT_TRANSFER_REQUESTED_M] = transfer_requested_m;
            record->limiter_inputs[INPUT_THETA_UPPER]          = theta_upper;
            record->limiter_inputs[INPUT_THETA_LOWER]          = theta_lower;
            record->limiter_inputs[INPUT_DZ_UPPER_M]           = dz_upper_m;
            record->limiter_inputs[INPUT_DZ_LOWER_M]           = dz_lower_m;
            record->limiter_inputs[INPUT_THETA_SAT]            = theta_sat;
            record->limiter_inputs[INPUT_PHI_SAT_CM]           = phi_sat_cm;
            record->limiter_inputs[INPUT_B_EXP]                = b_exp;
            record->limited_transfer_m = limited_transfer_m;
            n_activations++;
        } else {
            n_activations_dropped++;
        }
    }
    return limited_transfer_m;
}


/* ---- forcing ---- */
static double rain_mm_per_h[MAX_HOURS];
static double pet_mm_per_h[MAX_HOURS];

static int read_forcing_csv(const char *path, long *n_hours_out)
{
    FILE *fp;
    char line[1024];
    long n_hours = 0;
    int year, month, day, hour, minute, second;
    double rain_value, pet_value;

    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open forcing file %s\n", path);
        return 1;
    }
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "ERROR: forcing file %s is empty\n", path);
        fclose(fp);
        return 1;
    }
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        if (sscanf(line, "%d/%d/%d %d:%d:%d,%lf,%lf",
                   &year, &month, &day, &hour, &minute, &second,
                   &rain_value, &pet_value) != 8) {
            fprintf(stderr, "ERROR: malformed forcing record: %s", line);
            fclose(fp);
            return 1;
        }
        if (n_hours >= MAX_HOURS) {
            fprintf(stderr, "ERROR: more than %d forcing records\n", MAX_HOURS);
            fclose(fp);
            return 1;
        }
        rain_mm_per_h[n_hours] = rain_value;
        pet_mm_per_h[n_hours] = pet_value;
        n_hours++;
    }
    fclose(fp);
    *n_hours_out = n_hours;
    return 0;
}


/* ---- soil properties derived as in the driver and E10 ---- */
static double phi_sat_cm_from_ksat(double K_sat_cm_per_h)
{
    return PHI_SAT_FROM_KSAT_COEFFICIENT_CM * pow(K_sat_cm_per_h, PHI_SAT_FROM_KSAT_EXPONENT);
}

static double theta_at_pressure_ratio(double pressure_ratio, double theta_sat,
                                      double phi_sat_cm, double b_exp)
{
    double head_cm = pressure_ratio * ATMOSPHERIC_PRESSURE_HEAD_CM;
    if (head_cm <= phi_sat_cm) return theta_sat;
    return theta_sat * pow(phi_sat_cm / head_cm, 1.0 / b_exp);
}


/* ---- step 1: run the record and record the activations ---- */
static int run_record(long n_hours, double klf_m_per_h, double perc_limiter,
                      double K_sat_cm_per_h, double b_exp,
                      double *volume_balance_residual_m)
{
    SoilControl ctrl;
    SoilGeometry geom;
    SoilParameters params;
    SoilStateIn sin;
    SoilStateOut sout;
    SoilForcing forcing;
    SoilFluxes flux;
    TimestepSoilVolumeBalance volbal;
    double theta_current[NDISC];
    int ch_lut_hint[NDISC];
    double water_table_depth_m;
    double residual_sum_m = 0.0;
    long i_hour;
    int i_disc;

    memset(&ctrl, 0, sizeof(ctrl));
    ctrl.ndisc = NDISC;
    ctrl.deepest_root_disc = NDISC;
    ctrl.use_ch_lookup_table = 0;
    ctrl.apply_fc_perc_threshold = 0;
    ctrl.n_sub_setting = N_SUB_SETTING_DEFAULT;   /* fixed 4 substeps */
    ctrl.lateral_scheme = LATERAL_SCHEME_EXPONENTIAL;
    ctrl.rain_input_scheme = RAIN_INPUT_SCHEME_CASCADE;  /* infiltration flux upper boundary, TRUE (FLO 11/06/2026) */
    ctrl.noahmp_fidelity = NOAHMP_FIDELITY_AS_WRITTEN;  /* unused by DSBM */
    ctrl.dt_hours = 1.0;

    geom.dz_m[0] = 0.10;  geom.dz_m[1] = 0.30;  geom.dz_m[2] = 0.60;  geom.dz_m[3] = 1.00;
    geom.zc_m[0] = 0.05;  geom.zc_m[1] = 0.25;  geom.zc_m[2] = 0.70;  geom.zc_m[3] = 1.50;

    memset(&params, 0, sizeof(params));
    params.theta_r = 0.0;
    params.theta_sat = SOIL_THETA_SAT;
    params.K_sat_cm_per_h = K_sat_cm_per_h;
    params.phi_sat_cm = phi_sat_cm_from_ksat(K_sat_cm_per_h);
    params.b_exp = b_exp;
    params.perc_limiter_0_to_1 = perc_limiter;
    params.klf_m_per_h = klf_m_per_h;
    params.theta_fc = theta_at_pressure_ratio(FIELD_CAPACITY_PRESSURE_RATIO,
                                           params.theta_sat, params.phi_sat_cm, params.b_exp);
    params.theta_aet_eq_pet = params.theta_fc;
    params.theta_wp = theta_at_pressure_ratio(WILTING_POINT_PRESSURE_RATIO,
                                           params.theta_sat, params.phi_sat_cm, params.b_exp);
    // derived wilting point clamped at the dry-end limit, as in the driver
    if (params.theta_wp < CH_SCALED_SOIL_MOISTURE_MIN * params.theta_sat) {
        params.theta_wp = CH_SCALED_SOIL_MOISTURE_MIN * params.theta_sat;
    }

    printf("soil: theta_sat %.4f  phi_sat %.4f cm  theta_fc %.5f  theta_wp %.5f\n",
           params.theta_sat, params.phi_sat_cm, params.theta_fc, params.theta_wp);

    initialize_hydrostatic_from_storage(SOIL_DEPTH_M, params.theta_sat, params.phi_sat_cm,
                                        params.b_exp, geom.zc_m, INITIAL_STORAGE_M,
                                        &water_table_depth_m, theta_current);
    for (i_disc = 0; i_disc < NDISC; i_disc++) ch_lut_hint[i_disc] = -1;

    for (i_hour = 0; i_hour < n_hours; i_hour++) {
        current_hour_index = i_hour;
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            sin.theta_in[i_disc] = theta_current[i_disc];
            sin.ch_lut_hint_in[i_disc] = ch_lut_hint[i_disc];
        }
        compute_props_with_option_stateless(NULL, theta_current, sin.psi_in, sin.K_in,
                                            params.theta_r, params.theta_sat, params.K_sat_cm_per_h,
                                            params.phi_sat_cm, params.b_exp, ch_lut_hint);
        forcing.rain_mm_per_h = rain_mm_per_h[i_hour];
        forcing.pet_mm_per_h = pet_mm_per_h[i_hour];

        if (soil_step_one_hour_stateless(&ctrl, &geom, &params, NULL, &sin, &forcing,
                                         &sout, &flux, &volbal) != 0) {
            fprintf(stderr, "ERROR: DSBM kernel failed at hour %ld\n", i_hour);
            return 1;
        }
        residual_sum_m += volbal.residual_m;
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            theta_current[i_disc] = sout.theta_out[i_disc];
            ch_lut_hint[i_disc] = sout.ch_lut_hint_out[i_disc];
        }
    }
    *volume_balance_residual_m = residual_sum_m;
    return 0;
}


/* ---- step 2: derivatives at one activation ---- */

/* interface number (1..NDISC-1, between disc k and k+1), from the thicknesses */
static int interface_from_thicknesses(double dz_upper_m, double dz_lower_m)
{
    const double dz_m[NDISC] = {0.10, 0.30, 0.60, 1.00};
    int k;
    for (k = 0; k < NDISC - 1; k++) {
        if (fabs(dz_upper_m - dz_m[k]) < 1.0e-12 && fabs(dz_lower_m - dz_m[k + 1]) < 1.0e-12) {
            return k + 1;
        }
    }
    return 0;
}

static double limiter_primal(const double limiter_inputs[N_LIMITER_INPUTS])
{
    return limit_downward_transfer_to_equal_head(
        limiter_inputs[INPUT_TRANSFER_REQUESTED_M],
        limiter_inputs[INPUT_THETA_UPPER], limiter_inputs[INPUT_THETA_LOWER],
        limiter_inputs[INPUT_DZ_UPPER_M], limiter_inputs[INPUT_DZ_LOWER_M],
        limiter_inputs[INPUT_THETA_SAT], limiter_inputs[INPUT_PHI_SAT_CM],
        limiter_inputs[INPUT_B_EXP]);
}

/* which entry of limiter_inputs[] each active input is */
static const int active_to_input[N_LIMITER_ACTIVE] = {
    INPUT_TRANSFER_REQUESTED_M, INPUT_THETA_UPPER, INPUT_THETA_LOWER,
    INPUT_THETA_SAT, INPUT_PHI_SAT_CM, INPUT_B_EXP
};

static const char *active_name[N_LIMITER_ACTIVE] = {
    "transfer_requested_m", "theta_upper", "theta_lower",
    "theta_sat", "phi_sat_cm", "b_exp"
};

/* derivative units: m of transfer per unit of each active input */
static const char *active_units[N_LIMITER_ACTIVE] = {
    "m m-1", "m (m3 m-3)-1", "m (m3 m-3)-1", "m (m3 m-3)-1", "m cm-1", "m"
};

/*
 * Analytic implicit-function derivatives at the limited transfer.
 * Returns 1 if theta_upper_after or theta_lower_after is at a clamp of
 * psi_from_theta(), in which case the formulas do not apply.
 */
static int analytic_derivatives(const double limiter_inputs[N_LIMITER_INPUTS],
                                double limited_transfer_m,
                                double derivative_analytic[N_LIMITER_ACTIVE],
                                double *head_residual)
{
    double theta_upper = limiter_inputs[INPUT_THETA_UPPER];
    double theta_lower = limiter_inputs[INPUT_THETA_LOWER];
    double dz_upper_m = limiter_inputs[INPUT_DZ_UPPER_M];
    double dz_lower_m = limiter_inputs[INPUT_DZ_LOWER_M];
    double theta_sat = limiter_inputs[INPUT_THETA_SAT];
    double phi_sat_cm = limiter_inputs[INPUT_PHI_SAT_CM];
    double b_exp = limiter_inputs[INPUT_B_EXP];
    double distance_between_centers_m = 0.5 * (dz_upper_m + dz_lower_m);
    double theta_upper_after = theta_upper - limited_transfer_m / dz_upper_m;
    double theta_lower_after = theta_lower + limited_transfer_m / dz_lower_m;
    double psi_upper_m;
    double psi_lower_m;
    double lower_head_m;                 /* psi_lower + dz_int (m) */
    double dh_dV;                        /* (m-1) */
    double dh_dtheta_upper, dh_dtheta_lower, dh_dtheta_sat, dh_dphi_sat_cm, dh_db;
    int k;

    psi_upper_m = psi_from_theta(theta_upper_after, theta_sat, phi_sat_cm, b_exp);
    psi_lower_m = psi_from_theta(theta_lower_after, theta_sat, phi_sat_cm, b_exp);
    lower_head_m = psi_lower_m + distance_between_centers_m;
    *head_residual = log(psi_upper_m) - log(lower_head_m);

    for (k = 0; k < N_LIMITER_ACTIVE; k++) derivative_analytic[k] = 0.0;
    if (!(theta_upper_after > CH_SCALED_SOIL_MOISTURE_MIN * theta_sat && theta_upper_after < theta_sat &&
          theta_lower_after > CH_SCALED_SOIL_MOISTURE_MIN * theta_sat && theta_lower_after < theta_sat)) {
        return 1;
    }

    /* h = ln psi_upper - ln(psi_lower + dz_int), psi = phi (theta/theta_sat)^(-b) */
    dh_dV = b_exp / (theta_upper_after * dz_upper_m)
          + b_exp * psi_lower_m / (theta_lower_after * dz_lower_m * lower_head_m);
    dh_dtheta_upper = -b_exp / theta_upper_after;
    dh_dtheta_lower = b_exp * psi_lower_m / (theta_lower_after * lower_head_m);
    dh_dtheta_sat = (b_exp / theta_sat) * distance_between_centers_m / lower_head_m;
    dh_dphi_sat_cm = distance_between_centers_m / (phi_sat_cm * lower_head_m);
    dh_db = -log(theta_upper_after / theta_sat)
          + psi_lower_m * log(theta_lower_after / theta_sat) / lower_head_m;

    derivative_analytic[ACTIVE_TRANSFER_REQUESTED_M] = 0.0;
    derivative_analytic[ACTIVE_THETA_UPPER] = -dh_dtheta_upper / dh_dV;
    derivative_analytic[ACTIVE_THETA_LOWER] = -dh_dtheta_lower / dh_dV;
    derivative_analytic[ACTIVE_THETA_SAT] = -dh_dtheta_sat / dh_dV;
    derivative_analytic[ACTIVE_PHI_SAT_CM] = -dh_dphi_sat_cm / dh_dV;
    derivative_analytic[ACTIVE_B_EXP] = -dh_db / dh_dV;
    return 0;
}

static void tapenade_derivatives(const double limiter_inputs[N_LIMITER_INPUTS],
                                 double derivative_ad[N_LIMITER_ACTIVE],
                                 double *limited_transfer_ad_m)
{
    double seed[N_LIMITER_ACTIVE];
    int k, j;
    for (k = 0; k < N_LIMITER_ACTIVE; k++) {
        for (j = 0; j < N_LIMITER_ACTIVE; j++) seed[j] = 0.0;
        seed[k] = 1.0;
        equal_head_tangent_from_tapenade(limiter_inputs, seed,
                                         limited_transfer_ad_m, &derivative_ad[k]);
    }
}

static void finite_difference_derivatives(const double limiter_inputs[N_LIMITER_INPUTS],
                                          double derivative_fd[N_LIMITER_ACTIVE])
{
    double perturbed[N_LIMITER_INPUTS];
    double step, value_plus, value_minus;
    int k, j, input_index;
    for (k = 0; k < N_LIMITER_ACTIVE; k++) {
        input_index = active_to_input[k];
        for (j = 0; j < N_LIMITER_INPUTS; j++) perturbed[j] = limiter_inputs[j];
        step = FD_RELATIVE_STEP * fabs(limiter_inputs[input_index]);
        if (step <= 0.0) step = FD_RELATIVE_STEP;
        perturbed[input_index] = limiter_inputs[input_index] + step;
        value_plus = limiter_primal(perturbed);
        perturbed[input_index] = limiter_inputs[input_index] - step;
        value_minus = limiter_primal(perturbed);
        derivative_fd[k] = (value_plus - value_minus) / (2.0 * step);
    }
}

static double relative_error(double value, double reference)
{
    double scale = fabs(reference);
    if (scale < AD_ACCEPT_ABSOLUTE_DERIVATIVE) scale = AD_ACCEPT_ABSOLUTE_DERIVATIVE;
    return fabs(value - reference) / scale;
}


int main(int argc, char **argv)
{
    const char *output_dir;
    const char *forcing_path;
    double klf_m_per_h, perc_limiter, K_sat_cm_per_h, b_exp;
    long n_hours;
    double volume_balance_residual_m;
    char csv_path[2048];
    FILE *fp_csv;

    long n_newton_path = 0, n_zero_transfer = 0, n_at_clamp = 0;
    long n_newton_fail = 0, n_zero_fail = 0;
    long n_primal_mismatch = 0;     /* tangent's primal value differs from the production limiter */
    double max_rel_ad_vs_analytic[N_LIMITER_ACTIVE];
    double max_rel_fd_vs_analytic[N_LIMITER_ACTIVE];
    double max_abs_head_residual = 0.0;
    long worst_activation = -1;
    double worst_rel_error = 0.0;
    long i_act;
    int k;

    if (argc != 7) {
        fprintf(stderr, "Usage: %s output_directory forcing_csv_file k_lf_m_per_h perc_limiter K_sat_cm_per_h b_exp\n", argv[0]);
        return 1;
    }
    output_dir = argv[1];
    forcing_path = argv[2];
    klf_m_per_h = atof(argv[3]);
    perc_limiter = atof(argv[4]);
    K_sat_cm_per_h = atof(argv[5]);
    b_exp = atof(argv[6]);

    activations = (LimiterActivation *)malloc((size_t)MAX_ACTIVATIONS * sizeof(LimiterActivation));
    if (!activations) {
        fprintf(stderr, "ERROR: cannot allocate the activation records\n");
        return 1;
    }
    if (read_forcing_csv(forcing_path, &n_hours)) return 1;

    printf("E11 equal-head limiter derivative check\n");
    printf("parameters: k_lf %g m/h  perc_limiter %g  Ksat %g cm/h  b %g\n",
           klf_m_per_h, perc_limiter, K_sat_cm_per_h, b_exp);
    printf("forcing: %s, %ld hours\n", forcing_path, n_hours);

    if (run_record(n_hours, klf_m_per_h, perc_limiter, K_sat_cm_per_h, b_exp,
                   &volume_balance_residual_m)) return 1;
    printf("sum of hourly volume-balance residuals = %.4e m\n", volume_balance_residual_m);
    printf("limiter calls %ld, activations recorded %ld, dropped %ld\n",
           n_limiter_calls, n_activations, n_activations_dropped);

    snprintf(csv_path, sizeof(csv_path), "%s/activations.csv", output_dir);
    fp_csv = fopen(csv_path, "w");
    if (!fp_csv) {
        fprintf(stderr, "ERROR: cannot open %s\n", csv_path);
        return 1;
    }
    fprintf(fp_csv, "# hour_index,interface,transfer_requested_m,limited_transfer_m,head_residual_ln");
    for (k = 0; k < N_LIMITER_ACTIVE; k++) {
        fprintf(fp_csv, ",dV_d_%s_AD,dV_d_%s_analytic,dV_d_%s_FD",
                active_name[k], active_name[k], active_name[k]);
    }
    fprintf(fp_csv, "\n");

    for (k = 0; k < N_LIMITER_ACTIVE; k++) {
        max_rel_ad_vs_analytic[k] = 0.0;
        max_rel_fd_vs_analytic[k] = 0.0;
    }

    for (i_act = 0; i_act < n_activations; i_act++) {
        const LimiterActivation *record = &activations[i_act];
        double derivative_ad[N_LIMITER_ACTIVE];
        double derivative_analytic[N_LIMITER_ACTIVE];
        double derivative_fd[N_LIMITER_ACTIVE];
        double limited_transfer_ad_m;
        double head_residual;
        double error_here;
        int at_clamp;
        int failed = 0;

        tapenade_derivatives(record->limiter_inputs, derivative_ad, &limited_transfer_ad_m);
        if (fabs(limited_transfer_ad_m - record->limited_transfer_m) >
            1.0e-14 * fabs(record->limited_transfer_m)) {
            n_primal_mismatch++;
        }

        if (record->limited_transfer_m <= 0.0) {
            /* zero transfer: the exact derivative is zero for every input */
            n_zero_transfer++;
            for (k = 0; k < N_LIMITER_ACTIVE; k++) {
                if (derivative_ad[k] != 0.0) failed = 1;
            }
            if (failed) n_zero_fail++;
            continue;
        }

        n_newton_path++;
        at_clamp = analytic_derivatives(record->limiter_inputs, record->limited_transfer_m,
                                        derivative_analytic, &head_residual);
        if (at_clamp) {
            n_at_clamp++;
            continue;
        }
        finite_difference_derivatives(record->limiter_inputs, derivative_fd);
        if (fabs(head_residual) > max_abs_head_residual) max_abs_head_residual = fabs(head_residual);

        for (k = 0; k < N_LIMITER_ACTIVE; k++) {
            if (k == ACTIVE_TRANSFER_REQUESTED_M) {
                error_here = fabs(derivative_ad[k]);    /* exact value is 0 (m m-1) */
            } else {
                error_here = relative_error(derivative_ad[k], derivative_analytic[k]);
            }
            if (error_here > max_rel_ad_vs_analytic[k]) max_rel_ad_vs_analytic[k] = error_here;
            if (error_here > AD_ACCEPT_RELATIVE_ERROR) failed = 1;
            if (error_here > worst_rel_error) {
                worst_rel_error = error_here;
                worst_activation = i_act;
            }
            if (k != ACTIVE_TRANSFER_REQUESTED_M) {
                error_here = relative_error(derivative_fd[k], derivative_analytic[k]);
                if (error_here > max_rel_fd_vs_analytic[k]) max_rel_fd_vs_analytic[k] = error_here;
            }
        }
        if (failed) n_newton_fail++;

        fprintf(fp_csv, "%ld,%d,%.17g,%.17g,%.6e", record->hour_index,
                interface_from_thicknesses(record->limiter_inputs[INPUT_DZ_UPPER_M],
                                           record->limiter_inputs[INPUT_DZ_LOWER_M]),
                record->limiter_inputs[INPUT_TRANSFER_REQUESTED_M],
                record->limited_transfer_m, head_residual);
        for (k = 0; k < N_LIMITER_ACTIVE; k++) {
            fprintf(fp_csv, ",%.12e,%.12e,%.12e",
                    derivative_ad[k], derivative_analytic[k], derivative_fd[k]);
        }
        fprintf(fp_csv, "\n");
    }
    fclose(fp_csv);

    printf("\nactivations: Newton path %ld (of which %ld at a psi clamp, not checked), zero transfer %ld\n",
           n_newton_path, n_at_clamp, n_zero_transfer);
    printf("largest |equal-head residual| at the limited transfer = %.3e (ln head)\n",
           max_abs_head_residual);
    printf("\nNewton path: largest error of AD and of FD against the analytic\n");
    printf("implicit-function derivative (relative; for transfer_requested_m,\n");
    printf("absolute, since the exact derivative is 0)\n");
    printf("   %-22s %-14s %12s %12s\n", "input", "units", "AD", "FD");
    for (k = 0; k < N_LIMITER_ACTIVE; k++) {
        if (k == ACTIVE_TRANSFER_REQUESTED_M) {
            printf("   %-22s %-14s %12.3e %12s\n", active_name[k], active_units[k],
                   max_rel_ad_vs_analytic[k], "-");
        } else {
            printf("   %-22s %-14s %12.3e %12.3e\n", active_name[k], active_units[k],
                   max_rel_ad_vs_analytic[k], max_rel_fd_vs_analytic[k]);
        }
    }
    if (worst_activation >= 0) {
        printf("worst AD error %.3e at activation %ld (hour index %ld)\n",
               worst_rel_error, worst_activation, activations[worst_activation].hour_index);
    }
    printf("\nNewton-path activations failing (AD error > %.0e): %ld of %ld checked\n",
           AD_ACCEPT_RELATIVE_ERROR, n_newton_fail, n_newton_path - n_at_clamp);
    printf("zero-transfer activations with nonzero AD derivative: %ld of %ld\n",
           n_zero_fail, n_zero_transfer);
    printf("activations where the tangent's primal value differs from the production\n");
    printf("limiter by more than 1e-14 relative: %ld\n", n_primal_mismatch);
    printf("details: %s\n", csv_path);

    free(activations);
    if (n_newton_fail != 0 || n_zero_fail != 0 || n_activations_dropped != 0 ||
        n_primal_mismatch != 0) {
        printf("RESULT: FAIL\n");
        return 2;
    }
    printf("RESULT: PASS\n");
    return 0;
}
