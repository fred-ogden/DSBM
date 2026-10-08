/*
 * tapenade_soil_cost_primal.c
 *
 * Primal (undifferentiated) cost-function wrapper for the E10 Tapenade
 * experiments: four calibration parameters, two of them soil hydraulic
 * properties.
 *
 * Head function handed to Tapenade:
 *
 *   independents:  klf_m_per_h  perc_limiter_0_to_1  K_sat_cm_per_h  b_exp
 *   dependent:     cost_function_value
 *
 * phi_sat is NOT a calibration parameter.  It is computed from Ksat inside
 * the differentiated code with the same regression CFE3.1 uses when its
 * config file has soil_saturated_capillary_head_calc_from_ksat=TRUE (and
 * the DSBM driver reads the same keyword):
 *
 *   phi_sat (cm) = 10.415 * Ksat (cm/h) ^ (-0.3266)
 *
 * Ksat and phi_sat are inversely related soil properties; tying phi_sat to
 * Ksat keeps a calibration from proposing large values of both (physically
 * implausible, and a cause of instability seen with DDS in CFE), and it
 * removes one calibration parameter.
 *
 * theta_fc is then computed from phi_sat and b with the Clapp-Hornberger
 * retention curve at the field-capacity pressure head
 * (soil_field_capacity_Pcap_over_Patm_0_1 times atmospheric pressure head),
 * and theta_aet_eq_pet = theta_fc, as in src/bmi_soil_driver.c.  The
 * wilting point theta_wp is computed the same way at a capillary pressure
 * of 15 atmospheres (CFE3.1 WILTING_POINT_PCAP_OVER_PATM; the driver's
 * default).  Because all of these are computed inside the head function,
 * Tapenade carries the whole chain to the gradient:
 *
 *   Ksat -> K(theta)                                    (direct)
 *   Ksat -> phi_sat -> psi(theta), theta_fc, theta_wp   (through the regression)
 *   b    -> psi(theta), K(theta), theta_fc, theta_wp
 *
 * The cost function is the E8 one: SSE_lateral/SST_lateral +
 * SSE_percolation/SST_percolation (sum of 1 - NSE for hourly lateral flow
 * and percolation).  This is a twin experiment that tests the derivative
 * machinery only: the "observations" are DSBM's own outputs at known true
 * parameter values.  It says nothing about whether these parameters are
 * identifiable from discharge at a gauge, which needs CFE3.1 plus routing.
 *
 * Design rules carried over from E8: active parameters are plain scalar
 * arguments; no structure crosses the wrapper boundary; the time loop is
 * inside the wrapper; straight-line code inside the time loop; all locals
 * declared at the top of each function and no ternary operators
 * (required by Tapenade reverse mode).
 *
 * Other soil parameters are the configs/soil_params.dat values.  Analytic
 * Clapp-Hornberger, no lookup table, no field-capacity threshold on
 * percolation, exponential lateral scheme (production default).
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


#define SOIL_COST_THETA_SAT                 0.439
#define SOIL_COST_FC_PRESSURE_RATIO         0.333
#define SOIL_COST_WP_PRESSURE_RATIO         15.0    /* CFE3.1 WILTING_POINT_PCAP_OVER_PATM */
#define SOIL_COST_ATM_PRESSURE_HEAD_CM      1033.2274528

/* Ksat to phi_sat regression (CFE3.1 soil_saturated_capillary_head_calc_from_ksat) */
#define PHI_SAT_FROM_KSAT_COEFFICIENT_CM    10.415
#define PHI_SAT_FROM_KSAT_EXPONENT          (-0.3266)


/*
 * Saturated capillary head (cm) from saturated hydraulic conductivity
 * (cm/h).  Differentiated: dphi_sat/dKsat = -0.3266 phi_sat / Ksat.
 */
double soil_cost_phi_sat_cm_from_ksat(double K_sat_cm_per_h)
{
    double phi_sat_cm;

    phi_sat_cm = PHI_SAT_FROM_KSAT_COEFFICIENT_CM *
                 pow(K_sat_cm_per_h, PHI_SAT_FROM_KSAT_EXPONENT);
    return phi_sat_cm;
}


/*
 * Field-capacity soil moisture (m3/m3) from the Clapp-Hornberger retention
 * curve, theta = theta_sat (phi_sat / h)^(1/b), at h = field-capacity
 * pressure head.  If the field-capacity head is at or below phi_sat the
 * soil is saturated there and theta_fc = theta_sat (the same rule as
 * src/bmi_soil_driver.c).  With the regression above, phi_sat < 344 cm
 * (the field-capacity head here) for any Ksat above about 2e-5 cm/h, so
 * that branch does not arise within the calibration bounds.
 */
double soil_cost_theta_fc(double theta_sat, double phi_sat_cm, double b_exp)
{
    double field_capacity_head_cm;
    double theta_fc;

    field_capacity_head_cm = SOIL_COST_FC_PRESSURE_RATIO * SOIL_COST_ATM_PRESSURE_HEAD_CM;

    if (field_capacity_head_cm <= phi_sat_cm) {
        theta_fc = theta_sat;
    } else {
        theta_fc = theta_sat * pow(phi_sat_cm / field_capacity_head_cm, 1.0 / b_exp);
    }
    return theta_fc;
}


/*
 * Wilting-point soil moisture (m3/m3) from the Clapp-Hornberger retention
 * curve at a capillary pressure of 15 atmospheres, as in CFE3.1 and
 * src/bmi_soil_driver.c:  theta_wp = theta_sat (phi_sat / h_wp)^(1/b).
 * If the wilting-point head is at or below phi_sat, theta_wp = theta_sat
 * (does not arise for any realistic phi_sat: h_wp is about 15500 cm).
 */
double soil_cost_theta_wp(double theta_sat, double phi_sat_cm, double b_exp)
{
    double wilting_point_head_cm;
    double theta_wp;

    wilting_point_head_cm = SOIL_COST_WP_PRESSURE_RATIO * SOIL_COST_ATM_PRESSURE_HEAD_CM;

    if (wilting_point_head_cm <= phi_sat_cm) {
        theta_wp = theta_sat;
    } else {
        theta_wp = theta_sat * pow(phi_sat_cm / wilting_point_head_cm, 1.0 / b_exp);
    }
    return theta_wp;
}


/*
 * dsbm_soil_cost_from_params: the head function differentiated by Tapenade.
 *
 * Inputs:
 *   klf_m_per_h              ACTIVE: lateral subsurface flow rate constant (m/h)
 *   perc_limiter_0_to_1      ACTIVE: percolation limiter (-)
 *   K_sat_cm_per_h           ACTIVE: saturated hydraulic conductivity (cm/h)
 *   b_exp                    ACTIVE: Clapp-Hornberger pore-size exponent (-)
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
void dsbm_soil_cost_from_params(double klf_m_per_h,
                                double perc_limiter_0_to_1,
                                double K_sat_cm_per_h,
                                double b_exp,
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
    double lateral_error_m;
    double percolation_error_m;
    int kernel_status;
    int kernel_status_sum;
    int i_disc;
    int i_step;

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
    par.theta_sat = SOIL_COST_THETA_SAT;
    par.K_sat_cm_per_h = K_sat_cm_per_h;
    par.phi_sat_cm = soil_cost_phi_sat_cm_from_ksat(K_sat_cm_per_h);
    par.b_exp = b_exp;
    par.perc_limiter_0_to_1 = perc_limiter_0_to_1;
    par.klf_m_per_h = klf_m_per_h;
    par.theta_fc = soil_cost_theta_fc(par.theta_sat, par.phi_sat_cm, par.b_exp);
    par.theta_aet_eq_pet = par.theta_fc;
    par.theta_wp = soil_cost_theta_wp(par.theta_sat, par.phi_sat_cm, par.b_exp);

    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        theta_current[i_disc] = theta_in[i_disc];
        ch_lut_hint[i_disc] = -1;
    }
    sse_lateral_m2 = 0.0;
    sse_percolation_m2 = 0.0;
    kernel_status_sum = 0;

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
        kernel_status_sum = kernel_status_sum + kernel_status;

        lateral_this_hour_m = 0.0;
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            lateral_this_hour_m = lateral_this_hour_m + flux.lateral_by_disc_m[i_disc];
            theta_current[i_disc] = sout.theta_out[i_disc];
            ch_lut_hint[i_disc] = sout.ch_lut_hint_out[i_disc];
        }

        lateral_error_m = lateral_this_hour_m - observed_lateral_m[i_step];
        percolation_error_m = flux.percolation_to_gw_m - observed_percolation_m[i_step];
        sse_lateral_m2 = sse_lateral_m2 + lateral_error_m * lateral_error_m;
        sse_percolation_m2 = sse_percolation_m2 + percolation_error_m * percolation_error_m;
    }

    *cost_function_value = sse_lateral_m2 / sst_lateral_m2
                         + sse_percolation_m2 / sst_percolation_m2;
    if (kernel_status_sum != 0) {
        *cost_function_value = -1.0;
    }
}


/*
 * Not differentiated (not called by the head function): hourly lateral
 * flow and percolation series for a parameter set, used to make the
 * synthetic observations.  Same model setup as dsbm_soil_cost_from_params().
 */
void dsbm_soil_series_from_params(double klf_m_per_h,
                                  double perc_limiter_0_to_1,
                                  double K_sat_cm_per_h,
                                  double b_exp,
                                  const double *theta_in,
                                  int n_steps,
                                  int n_sub_setting,
                                  const double *rain_mm_per_h,
                                  const double *pet_mm_per_h,
                                  double *lateral_series_m,
                                  double *percolation_series_m)
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
    par.theta_sat = SOIL_COST_THETA_SAT;
    par.K_sat_cm_per_h = K_sat_cm_per_h;
    par.phi_sat_cm = soil_cost_phi_sat_cm_from_ksat(K_sat_cm_per_h);
    par.b_exp = b_exp;
    par.perc_limiter_0_to_1 = perc_limiter_0_to_1;
    par.klf_m_per_h = klf_m_per_h;
    par.theta_fc = soil_cost_theta_fc(par.theta_sat, par.phi_sat_cm, par.b_exp);
    par.theta_aet_eq_pet = par.theta_fc;
    par.theta_wp = soil_cost_theta_wp(par.theta_sat, par.phi_sat_cm, par.b_exp);

    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        theta_current[i_disc] = theta_in[i_disc];
        ch_lut_hint[i_disc] = -1;
    }

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

        soil_step_one_hour_stateless(
            &ctrl, &geom, &par, NULL, &sin, &forcing, &sout, &flux, &volbal);

        lateral_series_m[i_step] = 0.0;
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            lateral_series_m[i_step] = lateral_series_m[i_step] + flux.lateral_by_disc_m[i_disc];
            theta_current[i_disc] = sout.theta_out[i_disc];
            ch_lut_hint[i_disc] = sout.ch_lut_hint_out[i_disc];
        }
        percolation_series_m[i_step] = flux.percolation_to_gw_m;
    }
}
