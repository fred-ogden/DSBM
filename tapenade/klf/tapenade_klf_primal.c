/*
 * tapenade_klf_primal.c
 *
 * Primal (undifferentiated) wrapper for the k_lf Tapenade experiments.
 *
 * This is the ONLY file whose head function is handed to Tapenade:
 *
 *   independent:  klf_m_per_h
 *   dependents:   lateral_total_by_disc_m theta_out percolation_total_m
 *                 aet_total_m precipitation_excess_total_m
 *
 * Design rules that keep the generated tangent signature predictable:
 *
 *   1. The active calibration parameter is a plain scalar argument,
 *      klf_m_per_h, not a field of a caller-supplied SoilParameters.
 *
 *   2. No structure is passed in or out.  The entering soil moisture is
 *      a plain double array.  SoilStateIn, SoilParameters, etc. are all
 *      local to the wrapper, so any tangent structures Tapenade creates
 *      for them stay inside the generated code and never appear in the
 *      generated argument list.  (In the Ksat experiment a caller-supplied
 *      SoilStateIn produced an unexpected SoilStateIn_diff argument.)
 *
 *   3. The time loop is inside the wrapper.  n_steps = 1 is the
 *      single-timestep experiment; n_steps > 1 lets the static k_lf
 *      calibration parameter propagate through soil moisture and
 *      Darcy-Buckingham redistribution across timesteps.  One Tapenade
 *      run therefore serves every experiment.
 *
 * Between timesteps, psi and K are recomputed from the updated soil
 * moisture with the analytic Clapp-Hornberger relations, exactly as
 * src/bmi_soil_driver.c does.
 *
 * Soil parameters are the configs/soil_params.dat values, except that
 * k_lf is the argument.  Analytic Clapp-Hornberger, no lookup table,
 * no field-capacity threshold on percolation.  Adaptive n_sub is left
 * as in production unless n_sub_fixed > 0.
 *
 * Units:
 *   klf_m_per_h               m/h  (lateral rate constant, as in the kernel)
 *   theta                     m3/m3
 *   rain, PET                 mm/h
 *   lateral, percolation      m, summed over all n_steps
 *
 * ASCII only.
 */

#include <math.h>

#include "soil_data_types.h"
#include "soil_helpers.h"
#include "dsbm_soilmoisture_stateless.h"


/*
 * Baseline soil parameters (configs/soil_params.dat).
 * These are macros so the wrapper and the test driver cannot disagree.
 */
#define KLF_EXP_THETA_SAT                 0.439
#define KLF_EXP_K_SAT_CM_PER_H            1.2168
#define KLF_EXP_B_EXP                     4.05
#define KLF_EXP_PERC_LIMITER_0_TO_1       0.81
#define KLF_EXP_THETA_WP                  0.10
#define KLF_EXP_FC_PRESSURE_RATIO         0.333
#define KLF_EXP_ATM_PRESSURE_HEAD_CM      1033.2274528


/*
 * Saturated capillary head from Ksat, as in the production driver with
 * soil_saturated_capillary_head_calc_from_ksat=TRUE.  Ksat is not active
 * in this experiment, so this is a passive constant.
 */
double klf_experiment_phi_sat_cm(void)
{
    double phi_sat_cm;

    phi_sat_cm = 10.415 * pow(KLF_EXP_K_SAT_CM_PER_H, -0.3266);
    return phi_sat_cm;
}


/*
 * Field-capacity soil moisture (m3/m3) from the Clapp-Hornberger
 * retention curve at the field-capacity pressure head.  Same formula as
 * src/bmi_soil_driver.c, written without a ternary.
 */
double klf_experiment_theta_fc(void)
{
    double phi_sat_cm;
    double field_capacity_head_cm;
    double theta_fc;

    phi_sat_cm = klf_experiment_phi_sat_cm();
    field_capacity_head_cm =
        KLF_EXP_FC_PRESSURE_RATIO * KLF_EXP_ATM_PRESSURE_HEAD_CM;

    if (field_capacity_head_cm <= phi_sat_cm) {
        theta_fc = KLF_EXP_THETA_SAT;
    } else {
        theta_fc = KLF_EXP_THETA_SAT *
                   pow(phi_sat_cm / field_capacity_head_cm,
                       1.0 / KLF_EXP_B_EXP);
    }
    return theta_fc;
}


/*
 * dsbm_lateral_from_klf
 *
 * Inputs:
 *   klf_m_per_h            ACTIVE: lateral subsurface flow rate constant (m/h)
 *   theta_in[NDISC]        passive: entering soil moisture (m3/m3)
 *   n_steps                passive: number of hourly timesteps
 *   n_sub_fixed            passive: 0 = production adaptive substeps,
 *                          >0 = fixed substep count every timestep
 *   lateral_analytic       passive: 0 = forward-Euler lateral removal,
 *                          1 = exact exponential integration per substep
 *   n_sub_minimum          passive: floor on adaptive substeps (0 = none)
 *   lateral_substep_severity passive: 1 = adaptive substeps also limit the
 *                          lateral removal fraction per substep
 *   rain_mm_per_h[n_steps] passive: rainfall forcing (mm/h)
 *   pet_mm_per_h[n_steps]  passive: PET forcing (mm/h)
 *
 * Outputs:
 *   lateral_total_by_disc_m[NDISC]  DEPENDENT: lateral flow by disc,
 *                                   summed over all timesteps (m)
 *   theta_out[NDISC]                DEPENDENT: soil moisture after the
 *                                   last timestep (m3/m3)
 *   percolation_total_m             DEPENDENT: percolation to groundwater,
 *                                   summed over all timesteps (m)
 *   aet_total_m                     DEPENDENT: actual ET, all discs,
 *                                   summed over all timesteps (m)
 *   precipitation_excess_total_m    DEPENDENT: rain rejected at the soil
 *                                   surface, summed over timesteps (m)
 *
 *   The last two are included so the derivative itself can be checked
 *   for volume conservation:
 *     sum_i dz_i * d(theta_out_i)/dk
 *       + sum_i d(lateral_i)/dk + d(perc)/dk + d(AET)/dk + d(excess)/dk = 0
 *   n_sub_used_by_step[n_steps]     passive diagnostic: adaptive substep
 *                                   count used in each timestep
 */
void dsbm_lateral_from_klf(double klf_m_per_h,
                           const double *theta_in,
                           int n_steps,
                           int n_sub_fixed,
                           int lateral_analytic,
                           int n_sub_minimum,
                           int lateral_substep_severity,
                           const double *rain_mm_per_h,
                           const double *pet_mm_per_h,
                           double *lateral_total_by_disc_m,
                           double *theta_out,
                           double *percolation_total_m,
                           double *aet_total_m,
                           double *precipitation_excess_total_m,
                           int *n_sub_used_by_step)
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
    int i_disc;
    int i_step;
    int kernel_status;

    /* ---- control: analytic CH, no FC percolation threshold ---- */
    ctrl.ndisc = NDISC;
    ctrl.deepest_root_disc = NDISC;
    ctrl.use_ch_lookup_table = 0;
    ctrl.apply_fc_perc_threshold = 0;
    ctrl.dt_hours = 1.0;
    /* wrapper argument keeps its original meaning: 0 = adaptive, >0 = fixed */
    ctrl.n_sub_setting = n_sub_fixed;
    if (n_sub_fixed <= 0) {
        ctrl.n_sub_setting = N_SUB_SETTING_ADAPTIVE;
    }
    ctrl.n_sub_minimum = n_sub_minimum;
    ctrl.substep_lateral_severity = lateral_substep_severity;
    /* wrapper argument keeps its original meaning: 1 = exponential */
    if (lateral_analytic) {
        ctrl.lateral_scheme = LATERAL_SCHEME_EXPONENTIAL;
    } else {
        ctrl.lateral_scheme = LATERAL_SCHEME_FORWARD_EULER;
    }

    /* ---- Noah-MP disc geometry (m) ---- */
    geom.dz[0] = 0.10;
    geom.dz[1] = 0.30;
    geom.dz[2] = 0.60;
    geom.dz[3] = 1.00;

    geom.zc[0] = 0.05;
    geom.zc[1] = 0.25;
    geom.zc[2] = 0.70;
    geom.zc[3] = 1.50;

    /* ---- soil parameters; only klf_m_per_h is active ---- */
    par.theta_r = 0.0;
    par.theta_sat = KLF_EXP_THETA_SAT;
    par.K_sat_cm_per_h = KLF_EXP_K_SAT_CM_PER_H;
    par.phi_sat_cm = klf_experiment_phi_sat_cm();
    par.b_exp = KLF_EXP_B_EXP;
    par.perc_limiter_0_to_1 = KLF_EXP_PERC_LIMITER_0_TO_1;
    par.klf_m_per_h = klf_m_per_h;
    par.theta_fc = klf_experiment_theta_fc();
    par.theta_aet_eq_pet = par.theta_fc;
    par.theta_wp = KLF_EXP_THETA_WP;

    /* ---- initialize state and accumulators ---- */
    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        theta_current[i_disc] = theta_in[i_disc];
        ch_lut_hint[i_disc] = -1;
        lateral_total_by_disc_m[i_disc] = 0.0;
    }
    *percolation_total_m = 0.0;
    *aet_total_m = 0.0;
    *precipitation_excess_total_m = 0.0;

    /* ---- hourly time loop ---- */
    for (i_step = 0; i_step < n_steps; i_step++) {

        /* entering state for this timestep, psi and K from current theta */
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            sin.theta_in[i_disc] = theta_current[i_disc];
            sin.ch_lut_hint_in[i_disc] = ch_lut_hint[i_disc];
        }

        compute_props_with_option_stateless(
            NULL,
            theta_current,
            sin.psi_in,
            sin.K_in,
            par.theta_r,
            par.theta_sat,
            par.K_sat_cm_per_h,
            par.phi_sat_cm,
            par.b_exp,
            ch_lut_hint);

        forcing.rain_mm_per_h = rain_mm_per_h[i_step];
        forcing.pet_mm_per_h = pet_mm_per_h[i_step];

        kernel_status = soil_step_one_hour_stateless(
            &ctrl,
            &geom,
            &par,
            NULL,
            &sin,
            &forcing,
            &sout,
            &flux,
            &volbal);

        if (kernel_status != 0) {
            *percolation_total_m = -1.0;
            return;
        }

        /* accumulate step-integrated fluxes and advance the state */
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            lateral_total_by_disc_m[i_disc] =
                lateral_total_by_disc_m[i_disc] + flux.lateral_by_disc_m[i_disc];
            theta_current[i_disc] = sout.theta_out[i_disc];
            ch_lut_hint[i_disc] = sout.ch_lut_hint_out[i_disc];
        }
        *percolation_total_m = *percolation_total_m + flux.percolation_to_gw_m;
        *aet_total_m = *aet_total_m + volbal.AET_m;
        *precipitation_excess_total_m =
            *precipitation_excess_total_m + flux.surface_precipitation_excess_m;
        n_sub_used_by_step[i_step] = flux.n_sub_used;
    }

    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        theta_out[i_disc] = theta_current[i_disc];
    }
}
