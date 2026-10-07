/*
 * dsbm_soilmoisture_stateless.c
 *
 * Stateless Discrete Soil Moisture Balance Model (DSBM) calculation kernel.
 *
 * Discrete Soil Moisture Balance Model (DSBM)
 *
 * Author:
 *   Fred L. Ogden, Ph.D., P.E.
 *   NOAA/National Weather Service
 *
 * This software was developed by an employee of the United States
 * Government as part of official duties and is not subject to
 * copyright protection in the United States under 17 U.S.C. Section 105.
 *
 * License:
 *   Apache License, Version 2.0
 *   SPDX-License-Identifier: Apache-2.0
 *
 * See the repository LICENSE file for additional information.
 */

// ====================================================================
// Stateless discrete soil moisture balance module (DSBM) kernel) to 
// mimic the soil moisture flux routines use in Noah-MP.  Solves fluxes 
// between four soil discretizations (NDISC=4) that are the same thickness
// as those used in Noah-MP (0.1m, 0.3m, 0.6m, 1.0m).  My be used with
// other disc geometry, but won't be equivalent to Noah-MP.   
//
// This code calculates the Darcy-Buckingham flux between cells and applies
// a finite volume forward-Euler solution with time step controls.  This
// scheme is guaranteed to conserve mass.  It allows sinks due to plant
// water uptake and conceptual lateral flow.
//
// Positive vertical flux is downward.
// Units:
//   theta: m3/m3
//   psi:   m
//   K:     m/h
//   rainfall, PET inputs: mm/h (converted to m/h inside)
//   step-integrated totals: m
//   rates: m/h
// 
// Validated and found to mimic the soil moisture form of the Richards'
// equation used in Noah-MP with both NSE and KGE >0.998 on internal fluxes,
// soil moisture contents, and percolation flux out of the bottom of the
// soil domain.
//
// Author: Fred L. Ogden, Ph.D., P.E., Chief Scientist, NOAA/NWS OWP
//         National Water Center, Tuscaloosa, Alabama
//         2025-2026
// ====================================================================

#include <math.h>
#include "dsbm_soilmoisture_stateless.h"
#include "soil_helpers.h"

#ifndef THETA_MIN
#define THETA_MIN 1.0e-03
#endif

static inline double storage_sum_ndisc(const double *theta, const double *dz)
{
    double s = 0.0;
    for (int i = 0; i < NDISC; i++) s += theta[i] * dz[i];
    return s;
}

static inline int any_disc_above(const double *theta, double thresh, int ndisc)
{
    for (int i = 0; i < ndisc; i++) if (theta[i] > thresh) return 1;
    return 0;
}

// Fraction of a disc's soil moisture above theta_fc that the lateral
// linear reservoir would remove per substep (linearized, a = k dt/(dz
// (theta_sat - theta_fc))) at which the optional lateral severity equals 1.
// Matches the 0.20 used for the rain criterion below.
#define LATERAL_SEVERITY_FRACTION_PER_SUBSTEP 0.20

// Pick a conservative substep count based on initial fluxes and rainfall demand.
// Generic across NDISC.  lateral_severity is 0 unless the optional lateral
// criterion is enabled (ctrl->substep_lateral_severity).
static int choose_n_sub_generic(double dt_hours,
                                double rain_mm_per_h,
                                const double *dz,
                                const double *theta, double theta_sat,
                                const double *q0_m_per_h, // [0..NDISC-2]
                                int ndisc,
                                double lateral_severity)
{
    int nintf = ndisc - 1;

    // Max interface magnitude
    double qmax = 0.0;
    for (int i = 0; i < nintf; i++) {
        double a = fabs(q0_m_per_h[i]);
        if (a > qmax) qmax = a;
    }

    // Thinnest layer
    double dzmin = dz[0];
    for (int i = 1; i < ndisc; i++) if (dz[i] < dzmin) dzmin = dz[i];

    // Flux criterion: how much of the thinnest cell would move in dt_hours?
    double move_potential = qmax * dt_hours;                 // m
    double flux_ratio     = (dzmin > 0.0) ? (move_potential / (0.10 * dzmin)) : 0.0;

    // Rain criterion: fraction of top-cell storage capacity asked for this hour
    double rain_rate_m_per_h = rain_mm_per_h / 1000.0;
    double rain_hour_m       = rain_rate_m_per_h * dt_hours;
    double cap1_m            = (theta_sat - theta[0]) * dz[0];
    if (cap1_m < 1e-12) cap1_m = 1e-12;
    double rain_ratio        = rain_hour_m / (0.20 * cap1_m);

    double severity = 0.0;
    if (flux_ratio > severity) severity = flux_ratio;
    if (rain_ratio > severity) severity = rain_ratio;
    if (lateral_severity > severity) severity = lateral_severity;

    int n_sub;
    if (severity <= 1.0)       n_sub = (rain_mm_per_h > 0.0) ? 2 : 1;
    else if (severity <= 2.0)  n_sub = 4;
    else if (severity <= 3.0)  n_sub = 6;
    else if (severity <= 5.0)  n_sub = 8;
    else                       n_sub = 12;

    if (n_sub < 1)  n_sub = 1;
    if (n_sub > 12) n_sub = 12;
    return n_sub;
}

int soil_step_one_hour_stateless(
    const SoilControl        *ctrl,
    const SoilGeometry       *geom,
    const SoilParameters     *par,
    const SoilLookupTables   *lut,          // may be NULL
    const SoilStateIn        *sin,
    const SoilForcing        *forcing,
    SoilStateOut             *sout,
    SoilFluxes               *flux,
    TimestepSoilVolumeBalance      *volbal)
{

    // ---- local working state ------------------------------------------------
    double theta[NDISC];
    for (int i = 0; i < NDISC; i++) theta[i] = sin->theta_in[i];

    // outputs zeroed
    for (int i = 0; i < NDISC; i++) {
        flux->AET_by_disc_m[i] = 0.0;
        flux->lateral_by_disc_m[i] = 0.0;
        flux->interface_vol_m[i] = 0.0;
        flux->interface_rate_m_per_h[i] = 0.0;
        sout->ch_lut_hint_out[i] = sin->ch_lut_hint_in[i]; // start with input hints
    }
    flux->percolation_to_gw_m = 0.0;
    flux->rain_into_soil_m    = 0.0;
    flux->surface_precipitation_excess_m       = 0.0;
    flux->n_sub_used          = 0;

    volbal->in_rain_m       = 0.0;
    volbal->excess_m        = 0.0;
    volbal->perc_m          = 0.0;
    volbal->AET_m           = 0.0;
    volbal->lateral_m       = 0.0;
    volbal->delta_storage_m = 0.0;
    volbal->residual_m      = 0.0;

    const int ndisc  = ctrl->ndisc;
    const int nintf  = ndisc - 1;
    const double dtH = ctrl->dt_hours;

    const double rain_rate_m_per_h = forcing->rain_mm_per_h / 1000.0;
    const double pet_rate_m_per_h  = forcing->pet_mm_per_h  / 1000.0;

    double theta_floor;

    theta_floor = fmax(THETA_MIN, par->theta_r);

    // A) storage at start
    double storage_start;
    storage_start = storage_sum_ndisc(theta, geom->dz);

    // B) initial fluxes and n_sub
    double q0[NDISC > 1 ? NDISC-1 : 1];
    for (int i = 0; i < nintf; i++) {
        q0[i] = flux_DB_pair(sin->psi_in[i], sin->K_in[i],
                             sin->psi_in[i+1], sin->K_in[i+1],
                             geom->dz[i], geom->dz[i+1]);
    }
    // Substeps per timestep (ctrl->n_sub_setting):
    //   N_SUB_SETTING_DEFAULT (0): fixed at N_SUB_DEFAULT_COUNT (4).  A
    //     fixed count avoids the jumps in outputs and derivatives that the
    //     adaptive choice causes when a calibration parameter moves a
    //     timestep across a severity threshold (Tapenade experiments E5-E7).
    //   N_SUB_SETTING_ADAPTIVE (-1): the original adaptive choice below.
    //   > 0: fixed at that count.
    //
    // Optional refinements of the adaptive choice only (both off by default):
    //   ctrl->substep_lateral_severity: also require that the lateral
    //     linear reservoir removes no more than
    //     LATERAL_SEVERITY_FRACTION_PER_SUBSTEP of a disc's moisture above
    //     theta_fc per substep (linearized), for discs above theta_fc;
    //   ctrl->n_sub_minimum: a floor on the adaptive substep count.
    double lateral_severity = 0.0;
    if (ctrl->substep_lateral_severity && par->klf_m_per_h > 0.0) {
        double fc_to_sat_range = par->theta_sat - par->theta_fc;
        if (fc_to_sat_range < 1.0e-12) fc_to_sat_range = 1.0e-12;
        for (int i = 0; i < ndisc; i++) {
            if (theta[i] > par->theta_fc) {
                double lateral_exponent_per_hour =
                    par->klf_m_per_h * dtH / (geom->dz[i] * fc_to_sat_range);
                double disc_severity =
                    lateral_exponent_per_hour / LATERAL_SEVERITY_FRACTION_PER_SUBSTEP;
                if (disc_severity > lateral_severity) lateral_severity = disc_severity;
            }
        }
    }

    int n_sub;
    if (ctrl->n_sub_setting > 0) {
        n_sub = ctrl->n_sub_setting;
    } else if (ctrl->n_sub_setting == N_SUB_SETTING_DEFAULT) {
        n_sub = N_SUB_DEFAULT_COUNT;
    } else {
        n_sub = choose_n_sub_generic(dtH, forcing->rain_mm_per_h,
                                     geom->dz, theta, par->theta_sat, q0, ndisc,
                                     lateral_severity);
        if (n_sub < ctrl->n_sub_minimum) n_sub = ctrl->n_sub_minimum;
    }
    if (n_sub < 1) n_sub = 1;
    flux->n_sub_used = n_sub;

    const double dt_sub = dtH / (double)n_sub;

    // per-substep work arrays
    double psi[NDISC], K[NDISC];
    double store_cap[NDISC];
    double pot_downflux[NDISC > 1 ? NDISC-1 : 1];  // >=0
    double V_if[NDISC > 1 ? NDISC-1 : 1];          // signed desired
    double Accept[NDISC + 1];                      // [0..ndisc], nd = bottom
    for (int i = 0; i <= ndisc; i++) Accept[i] = 0.0;

    // External inflow (incident) is available to the caller via forcing and dtH.
    // Here we only track infiltrated vs excess for the step.
    // Loop over substeps
    for (int ss = 0; ss < n_sub; ss++) {
        // D1) rainfall into disc 0
        double rain_sub = rain_rate_m_per_h * dt_sub;  // m
        double cap1     = (par->theta_sat - theta[0]) * geom->dz[0];
        if (cap1 < 0.0) cap1 = 0.0;

        double used   = rain_sub;
        if (used > cap1) used = cap1;

        double excess = rain_sub - used;
        if (excess < 0.0) excess = 0.0;

        theta[0] += used / geom->dz[0];
        if (theta[0] > par->theta_sat) theta[0] = par->theta_sat;

        flux->rain_into_soil_m += used;
        flux->surface_precipitation_excess_m    += excess;

        // D2) properties after rainfall addition (LUT or analytic)
        int hint_local[NDISC];
        for (int i = 0; i < NDISC; i++) hint_local[i] = sout->ch_lut_hint_out[i];

        compute_props_with_option_stateless(
            (ctrl->use_ch_lookup_table ? lut : NULL),
            theta, psi, K,
            par->theta_r, par->theta_sat,
            par->K_sat_cm_per_h, par->phi_sat_cm, par->b_exp,
            hint_local);

        // keep updated hints
        for (int i = 0; i < NDISC; i++) sout->ch_lut_hint_out[i] = hint_local[i];

        // D3) interface fluxes and desired substep volumes
        for (int i = 0; i < nintf; i++) {
            double q;
            q = flux_DB_pair(psi[i], K[i], psi[i+1], K[i+1],
                                          geom->dz[i], geom->dz[i+1]);
            V_if[i] = q * dt_sub;
            pot_downflux[i] = (q > 0.0) ? (q * dt_sub) : 0.0;
        }

        // D4) per-disc free storage up to saturation
        for (int d = 0; d < ndisc; d++) {
            double s = (par->theta_sat - theta[d]) * geom->dz[d];
            if (s < 0.0) s = 0.0;
            store_cap[d] = s;
        }

        // D6) bottom potential percolation.
        // With --apply-fc-perc-threshold, preserve the original DSBM behavior:
        // drainage occurs only above theta_fc and cannot drain the bottom disc
        // below theta_fc.  Without it, use conductivity-based free drainage.
        double bottom_potential = 0.0;
        if (!ctrl->apply_fc_perc_threshold || theta[ndisc-1] > par->theta_fc) {
            double K_now = K_from_theta(theta[ndisc-1], par->theta_sat,
                                        par->K_sat_cm_per_h, par->b_exp);
            double percolation_rate = par->perc_limiter_0_to_1 * K_now; // m/h
            if (percolation_rate > 0.0) bottom_potential = percolation_rate * dt_sub;
        }

        // D7) downstream acceptance (bottom up)
        // Accept[ndisc] = last store + bottom_potential
        Accept[ndisc] = store_cap[ndisc-1] + bottom_potential;

        for (int i = nintf-1; i >= 0; i--) {
            double pass = pot_downflux[i];
            int down_index = i + 2;                 // downstream Accept slot; last interface -> ndisc
            if (down_index >= ndisc) down_index = ndisc;
            if (pass > Accept[down_index]) pass = Accept[down_index];
            Accept[i+1] = store_cap[i] + pass;
        }

        // D8) apply capped transfers across all interfaces
        for (int i = 0; i < nintf; i++) {
            double V = V_if[i];

            double accept_down =
                (i + 2 <= ndisc-1) ? Accept[i+2] : Accept[ndisc];

            if (V > 0.0) {
                double donor_avail = (theta[i] - theta_floor) * geom->dz[i];
                if (donor_avail < 0.0) donor_avail = 0.0;

                double recv_space = store_cap[i+1];

                double chain_pass = pot_downflux[i];
                if (chain_pass > accept_down) chain_pass = accept_down;

                double max_out = store_cap[i] + chain_pass;

                if (V > max_out) V = max_out;
                if (V > donor_avail) V = donor_avail;
                if (V > recv_space) V = recv_space; 
                if (V < 0.0) V = 0.0;

                theta[i]   -= V / geom->dz[i];
                theta[i+1] += V / geom->dz[i+1];

                if (theta[i]   < theta_floor)      theta[i]   = theta_floor;
                if (theta[i+1] > par->theta_sat)   theta[i+1] = par->theta_sat;
            }
            else if (V < 0.0) {
                double need = -V;

                double donor_avail = (theta[i+1] - theta_floor) * geom->dz[i+1];
                if (donor_avail < 0.0) donor_avail = 0.0;

                double recv_space = store_cap[i];

                double move = need;
                if (move > donor_avail) move = donor_avail;
                if (move > recv_space)  move = recv_space;

                V = -move;

                theta[i+1] -= move / geom->dz[i+1];
                theta[i]   += move / geom->dz[i];

                if (theta[i+1] < theta_floor)      theta[i+1] = theta_floor;
                if (theta[i]   > par->theta_sat)   theta[i]   = par->theta_sat;
            }

            // accumulate internal interface volume
            flux->interface_vol_m[i] += V;
        }

        // D9) apply bottom percolation
        double perc_vol = 0.0;
        if (bottom_potential > 0.0) {
            double drainage_floor = ctrl->apply_fc_perc_threshold
                                  ? par->theta_fc : theta_floor;
            double avail = (theta[ndisc-1] - drainage_floor) * geom->dz[ndisc-1];
            if (avail < 0.0) avail = 0.0;

            perc_vol = bottom_potential;
            if (perc_vol > avail) perc_vol = avail;

            theta[ndisc-1] -= perc_vol / geom->dz[ndisc-1];
            if (theta[ndisc-1] < drainage_floor) theta[ndisc-1] = drainage_floor;
        }
        flux->percolation_to_gw_m += perc_vol;
        flux->interface_vol_m[ndisc-1] += perc_vol;

        // D10) ET from root zone
        int nroot = ctrl->deepest_root_disc;
        if (nroot < 1) nroot = 1;
        if (nroot > ndisc) nroot = ndisc;

        double pet_sub = pet_rate_m_per_h * dt_sub;
        double denom = par->theta_aet_eq_pet - par->theta_wp;
        if (denom < 1.0e-12) denom = 1.0e-12;

        double et_this_sub = 0.0;

        /*
         * Distribute PET equally among the root-zone discs, consistent
         * with the Noah-MP soil-moisture formulation that DSBM was
         * developed to mimic.  Each disc independently realizes its
         * assigned fraction of PET according to its local soil-moisture
         * stress.
         *
         * Field capacity defines the moisture content at and above which
         * AET equals the assigned PET demand.  Between field capacity and
         * wilting point, AET decreases linearly to zero.  Using theta_fc
         * as the stress-onset threshold avoids an additional calibration
         * parameter.
         *
         * Extraction from any disc is limited so that theta cannot fall
         * below the wilting point during the substep.
         */
        for (int i = 0; i < nroot; i++) {
            double root_frac = 1.0 / (double)nroot;

            double th = theta[i];
            double f_aet;  // fraction of assigned PET realized
            if (th <= par->theta_wp)                  f_aet = 0.0;
            else if (th >= par->theta_aet_eq_pet)     f_aet = 1.0;
            else                                      f_aet = (th - par->theta_wp) / denom;

            double demand = f_aet * pet_sub * root_frac;

            double avail = (th - par->theta_wp) * geom->dz[i];
            if (avail < 0.0) avail = 0.0;

            double AET_i = demand;
            if (AET_i > avail) AET_i = avail;

            theta[i] -= AET_i / geom->dz[i];
            // Roundoff guard only.  AET_i <= avail keeps theta at or above
            // theta_wp when the disc started above it.  A disc that entered
            // ET already below theta_wp (drained there by Darcy-Buckingham
            // flow) has AET_i = 0 and must stay where it is: raising it to
            // theta_wp would create soil moisture with no matching flux.
            if (th > par->theta_wp && theta[i] < par->theta_wp) theta[i] = par->theta_wp;

            et_this_sub            += AET_i;
            flux->AET_by_disc_m[i] += AET_i;
        }

        // D11) lateral removal to Nash (only if any disc > FC)
        // Default: the linear-reservoir removal above theta_fc is integrated
        // exactly over the substep (exponential decay toward theta_fc; no
        // storage cap needed, no zero-gradient region for calibration).
        // ctrl->lateral_scheme = LATERAL_SCHEME_FORWARD_EULER restores the
        // legacy explicit step.
        if (par->klf_m_per_h > 0.0 && any_disc_above(theta, par->theta_fc, ndisc)) {
            double lat_removed;
            if (ctrl->lateral_scheme == LATERAL_SCHEME_FORWARD_EULER) {
                lat_removed =
                    remove_lateral_to_subsurface_nash_substep(
                        theta, geom->dz, par->theta_fc, par->theta_sat,
                        par->klf_m_per_h, dt_sub, flux->lateral_by_disc_m);
            } else {
                lat_removed =
                    remove_lateral_to_subsurface_nash_substep_exponential(
                        theta, geom->dz, par->theta_fc, par->theta_sat,
                        par->klf_m_per_h, dt_sub, flux->lateral_by_disc_m);
            }

            volbal->lateral_m += lat_removed;
        }

        // accumulate step totals
        volbal->AET_m += et_this_sub;

        // D12) safety clamp
        for (int i = 0; i < ndisc; i++) {
            if (theta[i] < theta_floor)    theta[i] = theta_floor;
            if (theta[i] > par->theta_sat) theta[i] = par->theta_sat;
        }
    } // end substeps

    // finalize rates
    for (int i = 0; i < NDISC; i++) {
        flux->interface_rate_m_per_h[i] = flux->interface_vol_m[i] / dtH;
    }

    // write state out
    for (int i = 0; i < NDISC; i++) sout->theta_out[i] = theta[i];

    // volume balance
    double storage_end;
    storage_end = storage_sum_ndisc(theta, geom->dz);

    volbal->in_rain_m  = flux->rain_into_soil_m;
    volbal->excess_m   = flux->surface_precipitation_excess_m;
    volbal->perc_m     = flux->percolation_to_gw_m;

    // total lateral already accumulated in volbal->lateral_m inside loop
    // Ensure it matches sum of per-disc laterals (defensive)
    double lat_sum = 0.0; for (int i = 0; i < NDISC; i++) lat_sum += flux->lateral_by_disc_m[i];
    volbal->lateral_m = lat_sum;

    volbal->delta_storage_m = (storage_end - storage_start);

    volbal->residual_m = (volbal->in_rain_m)
                   - (volbal->perc_m + volbal->AET_m + volbal->lateral_m)
                   -  volbal->delta_storage_m;

    return 0;
}
