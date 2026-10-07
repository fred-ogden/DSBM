/*
 * tapenade_cost_primal.c
 *
 * Primal (undifferentiated) cost-function wrapper for the E8 Tapenade
 * experiments: the first calibration gradient.
 *
 * Head function handed to Tapenade:
 *
 *   independents:  klf_m_per_h  perc_limiter_0_to_1
 *   dependent:     cost_function_value
 *
 * The cost function compares two simulated hourly series against
 * "observed" series of the same quantities:
 *
 *   lateral_m(t)     = lateral subsurface flow from all discs in hour t (m)
 *   percolation_m(t) = percolation to groundwater in hour t (m)
 *
 *   cost = SSE_lateral / SST_lateral  +  SSE_percolation / SST_percolation
 *
 * where SSE = sum over hours of (simulated - observed)^2 and SST = sum over
 * hours of (observed - mean observed)^2.  Each term is 1 - NSE, so the cost
 * is 0 for a perfect fit and the two series carry equal weight despite
 * their different magnitudes.  SST depends only on the observations, so it
 * is computed by the caller and passed in as a passive constant.
 *
 * Design rules carried over from the k_lf experiments:
 *   - active parameters are plain scalar arguments;
 *   - no structure crosses the wrapper boundary;
 *   - the time loop is inside the wrapper;
 *   - psi and K are recomputed from theta at the start of every hour,
 *     as in src/bmi_soil_driver.c.
 *
 * Soil parameters other than the two active ones are the
 * configs/soil_params.dat values.  Analytic Clapp-Hornberger, no lookup
 * table, no field-capacity threshold on percolation, exponential lateral
 * scheme (production default).
 *
 * n_sub_setting follows SoilControl.n_sub_setting:
 *   N_SUB_SETTING_DEFAULT (0) = fixed 4, N_SUB_SETTING_ADAPTIVE (-1),
 *   > 0 = fixed count.
 *
 * ASCII only.
 */

#include <math.h>

#include "soil_data_types.h"
#include "soil_helpers.h"
#include "dsbm_soilmoisture_stateless.h"


#define COST_EXP_THETA_SAT                 0.439
#define COST_EXP_K_SAT_CM_PER_H            1.2168
#define COST_EXP_B_EXP                     4.05
#define COST_EXP_THETA_WP                  0.10
#define COST_EXP_FC_PRESSURE_RATIO         0.333
#define COST_EXP_ATM_PRESSURE_HEAD_CM      1033.2274528


/* Saturated capillary head from Ksat (calc_from_ksat=TRUE); passive here. */
double cost_experiment_phi_sat_cm(void)
{
    double phi_sat_cm;

    phi_sat_cm = 10.415 * pow(COST_EXP_K_SAT_CM_PER_H, -0.3266);
    return phi_sat_cm;
}


/* Field-capacity soil moisture from the Clapp-Hornberger retention curve. */
double cost_experiment_theta_fc(void)
{
    double phi_sat_cm;
    double field_capacity_head_cm;
    double theta_fc;

    phi_sat_cm = cost_experiment_phi_sat_cm();
    field_capacity_head_cm =
        COST_EXP_FC_PRESSURE_RATIO * COST_EXP_ATM_PRESSURE_HEAD_CM;

    if (field_capacity_head_cm <= phi_sat_cm) {
        theta_fc = COST_EXP_THETA_SAT;
    } else {
        theta_fc = COST_EXP_THETA_SAT *
                   pow(phi_sat_cm / field_capacity_head_cm,
                       1.0 / COST_EXP_B_EXP);
    }
    return theta_fc;
}


/*
 * Run the soil column for n_steps hours.  Shared by the cost wrapper and
 * the series generator below.  If lateral_series_m / percolation_series_m
 * are non-NULL they receive the hourly series; if observed series are
 * non-NULL the cost is accumulated.
 */
static void dsbm_run_hours(double klf_m_per_h,
                           double perc_limiter_0_to_1,
                           const double *theta_in,
                           int n_steps,
                           int n_sub_setting,
                           const double *rain_mm_per_h,
                           const double *pet_mm_per_h,
                           const double *observed_lateral_m,
                           const double *observed_percolation_m,
                           double sst_lateral_m2,
                           double sst_percolation_m2,
                           double *lateral_series_m,
                           double *percolation_series_m,
                           double *cost_function_value)
{
    SoilControl ctrl;
    SoilGeometry geom;
    SoilParameters par;
    SoilStateIn sin;
    SoilStateOut sout;
    SoilForcing forcing;
    SoilFluxes flux;
    TimestepSoilVolumeBalance volbal;

    double theta_current[NDISC];
    int ch_lut_hint[NDISC];
    double sse_lateral_m2;
    double sse_percolation_m2;
    double lateral_this_hour_m;
    double error_m;
    int i_disc;
    int i_step;
    int kernel_status;

    ctrl.ndisc = NDISC;
    ctrl.deepest_root_disc = NDISC;
    ctrl.use_ch_lookup_table = 0;
    ctrl.apply_fc_perc_threshold = 0;
    ctrl.dt_hours = 1.0;
    ctrl.n_sub_setting = n_sub_setting;
    ctrl.n_sub_minimum = 0;
    ctrl.substep_lateral_severity = 0;
    ctrl.lateral_scheme = LATERAL_SCHEME_EXPONENTIAL;

    geom.dz[0] = 0.10;
    geom.dz[1] = 0.30;
    geom.dz[2] = 0.60;
    geom.dz[3] = 1.00;
    geom.zc[0] = 0.05;
    geom.zc[1] = 0.25;
    geom.zc[2] = 0.70;
    geom.zc[3] = 1.50;

    par.theta_r = 0.0;
    par.theta_sat = COST_EXP_THETA_SAT;
    par.K_sat_cm_per_h = COST_EXP_K_SAT_CM_PER_H;
    par.phi_sat_cm = cost_experiment_phi_sat_cm();
    par.b_exp = COST_EXP_B_EXP;
    par.perc_limiter_0_to_1 = perc_limiter_0_to_1;
    par.klf_m_per_h = klf_m_per_h;
    par.theta_fc = cost_experiment_theta_fc();
    par.theta_aet_eq_pet = par.theta_fc;
    par.theta_wp = COST_EXP_THETA_WP;

    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        theta_current[i_disc] = theta_in[i_disc];
        ch_lut_hint[i_disc] = -1;
    }
    sse_lateral_m2 = 0.0;
    sse_percolation_m2 = 0.0;

    for (i_step = 0; i_step < n_steps; i_step++) {
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            sin.theta_in[i_disc] = theta_current[i_disc];
            sin.ch_lut_hint_in[i_disc] = ch_lut_hint[i_disc];
        }
        compute_props_with_option_stateless(
            NULL, theta_current, sin.psi_in, sin.K_in,
            par.theta_r, par.theta_sat, par.K_sat_cm_per_h,
            par.phi_sat_cm, par.b_exp, ch_lut_hint);

        forcing.rain_mm_per_h = rain_mm_per_h[i_step];
        forcing.pet_mm_per_h = pet_mm_per_h[i_step];

        kernel_status = soil_step_one_hour_stateless(
            &ctrl, &geom, &par, NULL, &sin, &forcing, &sout, &flux, &volbal);
        if (kernel_status != 0) {
            *cost_function_value = -1.0;
            return;
        }

        lateral_this_hour_m = 0.0;
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            lateral_this_hour_m = lateral_this_hour_m + flux.lateral_by_disc_m[i_disc];
            theta_current[i_disc] = sout.theta_out[i_disc];
            ch_lut_hint[i_disc] = sout.ch_lut_hint_out[i_disc];
        }

        if (lateral_series_m != NULL) lateral_series_m[i_step] = lateral_this_hour_m;
        if (percolation_series_m != NULL) percolation_series_m[i_step] = flux.percolation_to_gw_m;

        if (observed_lateral_m != NULL) {
            error_m = lateral_this_hour_m - observed_lateral_m[i_step];
            sse_lateral_m2 = sse_lateral_m2 + error_m * error_m;
            error_m = flux.percolation_to_gw_m - observed_percolation_m[i_step];
            sse_percolation_m2 = sse_percolation_m2 + error_m * error_m;
        }
    }

    *cost_function_value = sse_lateral_m2 / sst_lateral_m2
                         + sse_percolation_m2 / sst_percolation_m2;
}


/*
 * dsbm_cost_from_params: the head function differentiated by Tapenade.
 *
 * Inputs:
 *   klf_m_per_h              ACTIVE: lateral subsurface flow rate constant (m/h)
 *   perc_limiter_0_to_1      ACTIVE: percolation limiter (-)
 *   theta_in[NDISC]          passive: entering soil moisture (m3/m3)
 *   n_steps                  passive: number of hourly timesteps
 *   n_sub_setting            passive: substeps (see header comment)
 *   rain_mm_per_h[n_steps]   passive forcing
 *   pet_mm_per_h[n_steps]    passive forcing
 *   observed_lateral_m[n_steps], observed_percolation_m[n_steps]
 *                            passive observations
 *   sst_lateral_m2, sst_percolation_m2
 *                            passive: sums of squared deviations of the
 *                            observations from their means
 * Output:
 *   cost_function_value      DEPENDENT: cost function (-)
 */
void dsbm_cost_from_params(double klf_m_per_h,
                           double perc_limiter_0_to_1,
                           const double *theta_in,
                           int n_steps,
                           int n_sub_setting,
                           const double *rain_mm_per_h,
                           const double *pet_mm_per_h,
                           const double *observed_lateral_m,
                           const double *observed_percolation_m,
                           double sst_lateral_m2,
                           double sst_percolation_m2,
                           double *cost_function_value)
{
    dsbm_run_hours(klf_m_per_h, perc_limiter_0_to_1, theta_in, n_steps,
                   n_sub_setting, rain_mm_per_h, pet_mm_per_h,
                   observed_lateral_m, observed_percolation_m,
                   sst_lateral_m2, sst_percolation_m2,
                   NULL, NULL, cost_function_value);
}


/*
 * Not differentiated: hourly series for a parameter set (used to make the
 * synthetic observations).
 */
void dsbm_series_from_params(double klf_m_per_h,
                             double perc_limiter_0_to_1,
                             const double *theta_in,
                             int n_steps,
                             int n_sub_setting,
                             const double *rain_mm_per_h,
                             const double *pet_mm_per_h,
                             double *lateral_series_m,
                             double *percolation_series_m)
{
    double unused_cost;

    dsbm_run_hours(klf_m_per_h, perc_limiter_0_to_1, theta_in, n_steps,
                   n_sub_setting, rain_mm_per_h, pet_mm_per_h,
                   NULL, NULL, 1.0, 1.0,
                   lateral_series_m, percolation_series_m, &unused_cost);
}
