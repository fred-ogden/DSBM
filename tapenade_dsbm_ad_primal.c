#include <stdio.h>
#include <math.h>

#include "soil_data_types.h"
#include "soil_helpers.h"
#include "dsbm_soilmoisture_stateless.h"

/*
 * One-timestep DSBM test for automatic differentiation.
 *
 * The entering soil state is initialized outside the differentiated
 * function and is therefore held fixed.  Within the differentiated
 * function, K_sat_cm_per_h is the only active input and
 * percolation_to_gw_m is the dependent output.
 */

void initialize_baseline_state_for_driver(SoilStateIn *sin)
{
    SoilGeometry geom;
    SoilParameters par;
    double theta[NDISC];
    double zwt;
    int hint[NDISC];

    geom.dz[0] = 0.10;
    geom.dz[1] = 0.30;
    geom.dz[2] = 0.60;
    geom.dz[3] = 1.00;

    geom.zc[0] = 0.05;
    geom.zc[1] = 0.25;
    geom.zc[2] = 0.70;
    geom.zc[3] = 1.50;

    par.theta_r = 0.0;
    par.theta_sat = 0.439;
    par.phi_sat_cm = 10.415 * pow(1.2168, -0.3266);
    par.b_exp = 4.05;

    for (int i = 0; i < NDISC; i++)
        hint[i] = -1;

    initialize_hydrostatic_from_storage(
        2.0,
        par.theta_sat,
        par.phi_sat_cm,
        par.b_exp,
        geom.zc,
        0.7843498367003,
        &zwt,
        theta);

    for (int i = 0; i < NDISC; i++) {
        sin->theta_in[i] = theta[i];
        sin->ch_lut_hint_in[i] = hint[i];
    }

    /*
     * Deliberately use the baseline Ksat here.  The complete entering
     * state, including K_in, is fixed while the timestep Ksat varies.
     */
    compute_props_with_option_stateless(
        NULL,
        theta,
        sin->psi_in,
        sin->K_in,
        par.theta_r,
        par.theta_sat,
        1.2168,
        par.phi_sat_cm,
        par.b_exp,
        hint);
}


void dsbm_percolation_from_ksat(double K_sat_cm_per_h,
                                const SoilStateIn *sin,
                                double *percolation_to_gw_m)
{
    SoilControl ctrl;
    SoilGeometry geom;
    SoilParameters par;
    SoilStateOut sout;
    SoilForcing forcing;
    SoilFluxes flux;
    TimestepSoilVolumeBalance volbal;

    ctrl.ndisc = NDISC;
    ctrl.deepest_root_disc = NDISC;
    ctrl.use_ch_lookup_table = 0;
    ctrl.apply_fc_perc_threshold = 0;
    ctrl.n_sub_fixed = 0;
    ctrl.lateral_scheme = LATERAL_SCHEME_FORWARD_EULER;  /* scheme used when this experiment was run */
    ctrl.dt_hours = 1.0;

    geom.dz[0] = 0.10;
    geom.dz[1] = 0.30;
    geom.dz[2] = 0.60;
    geom.dz[3] = 1.00;

    geom.zc[0] = 0.05;
    geom.zc[1] = 0.25;
    geom.zc[2] = 0.70;
    geom.zc[3] = 1.50;

    /*
     * Baseline soil parameters from configs/soil_params.dat.
     *
     * phi_sat_cm and theta_fc are deliberately held fixed in this first
     * experiment.  We are differentiating the DSBM timestep kernel with
     * respect to K_sat_cm_per_h, not the configuration preprocessing that
     * derives other parameters from Ksat.
     */
    par.theta_r = 0.0;
    par.theta_sat = 0.439;
    par.K_sat_cm_per_h = K_sat_cm_per_h;
    par.phi_sat_cm = 10.415 * pow(1.2168, -0.3266);
    par.b_exp = 4.05;
    par.perc_limiter_0_to_1 = 0.81;
    par.klf_m_per_h = 0.000005;

    {
        const double atmospheric_pressure_head_cm = 1033.2274528;
        const double field_capacity_pressure_ratio = 0.333;
        const double field_capacity_head_cm =
            field_capacity_pressure_ratio * atmospheric_pressure_head_cm;

        par.theta_fc =
            (field_capacity_head_cm <= par.phi_sat_cm)
            ? par.theta_sat
            : par.theta_sat *
              pow(par.phi_sat_cm / field_capacity_head_cm,
                  1.0 / par.b_exp);

        par.theta_aet_eq_pet = par.theta_fc;
    }

    par.theta_wp = 0.10;

    forcing.rain_mm_per_h = 0.0;
    forcing.pet_mm_per_h = 0.0;

    if (soil_step_one_hour_stateless(
            &ctrl,
            &geom,
            &par,
            NULL,
            sin,
            &forcing,
            &sout,
            &flux,
            &volbal) != 0) {
        *percolation_to_gw_m = -1.0;
        return;
    }

    *percolation_to_gw_m = flux.percolation_to_gw_m;
}
