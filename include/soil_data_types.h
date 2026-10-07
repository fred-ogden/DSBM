/*
 * soil_data_types.h
 *
 * Shared data structures for DSBM and Noah-MP soil-moisture calculations.
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

#ifndef SOIL_TYPES_H
#define SOIL_TYPES_H

#include "soil_config.h"

#ifdef __cplusplus
extern "C" {
#endif

// Lateral subsurface flow integration within each DSBM substep.
// The exponential value is 0 so that a zero-initialized SoilControl
// selects it.
#define LATERAL_SCHEME_EXPONENTIAL    0   // exact linear-reservoir solution
#define LATERAL_SCHEME_FORWARD_EULER  1   // legacy explicit step with storage cap

// Substeps per timestep.  The zero value of n_sub_setting selects the
// default fixed count, so a zero-initialized SoilControl gets it.
#define N_SUB_DEFAULT_COUNT       4   // default fixed substeps per timestep
#define N_SUB_SETTING_DEFAULT     0   // use N_SUB_DEFAULT_COUNT
#define N_SUB_SETTING_ADAPTIVE   (-1) // original adaptive substep choice

typedef struct {
    int    ndisc;                 // must equal NDISC
    int    deepest_root_disc;     // 1..ndisc
    int    use_ch_lookup_table;   // 1 => use LUT; 0 => analytic CH
    int    apply_fc_perc_threshold; // 1 => bottom drainage only above theta_fc
    int    n_sub_setting;         // N_SUB_SETTING_DEFAULT (0): fixed N_SUB_DEFAULT_COUNT;
                                  // N_SUB_SETTING_ADAPTIVE (-1): original adaptive choice;
                                  // >0: fixed at that many substeps per timestep
    int    n_sub_minimum;         // adaptive only: 0 => no floor; >0 => n_sub is at least this
    int    substep_lateral_severity; // adaptive only: 1 => also consider lateral removal rate
    int    lateral_scheme;        // LATERAL_SCHEME_EXPONENTIAL (0, default) or
                                  // LATERAL_SCHEME_FORWARD_EULER (1)
    double dt_hours;              // usually 1.0
} SoilControl;

typedef struct {
    double dz[NDISC];
    double zc[NDISC];
} SoilGeometry;

typedef struct {
    double theta_r;               // residual saturation (m3/m3)
    double theta_sat;             // saturation (m3/m3)
    double theta_fc;              // field capacity (m3/m3)
    double theta_wp;              // wilting point (m3/m3)
    double theta_aet_eq_pet;      // AET = PET at / above this \u03b8 (m3/m3)

    double K_sat_cm_per_h;        // cm/h
    double phi_sat_cm;            // cm
    double b_exp;                 // Clapp\u2013Hornberger exponent

    double perc_limiter_0_to_1;   // 0..1 bottom drainage limiter
    double klf_m_per_h;           // lateral removal rate constant (m/h)
} SoilParameters;

// CH lookup tables over Theta=(theta-theta_r)/(theta_sat-theta_r)
typedef struct {
    int    n;
    double lnTheta_min;
    double dlnTheta;
    double inv_dlnTheta;
    double *lnpsi;                // ln(psi[m]) length n
    double *lnK;                  // ln(K[m/h]) length n
} SoilLookupTables;

typedef struct {
     double theta_in[NDISC];
     double psi_in[NDISC];   // m
     double K_in[NDISC];     // m/h
     int    ch_lut_hint_in[NDISC];
} SoilStateIn;

typedef struct {
    double theta_out[NDISC];
    int    ch_lut_hint_out[NDISC];
} SoilStateOut;

typedef struct {
    double rain_mm_per_h;         // mm/h
    double pet_mm_per_h;          // mm/h
} SoilForcing;

typedef struct {
    // Step-integrated exchanges (m)
    double AET_by_disc_m[NDISC];
    double lateral_by_disc_m[NDISC];
    double percolation_to_gw_m;
    double rain_into_soil_m;
    double surface_precipitation_excess_m;

    // Internal vertical exchanges: [0..NDISC-2] interfaces i\u2192i+1; [NDISC-1] bottom perc
    double interface_vol_m[NDISC];
    double interface_rate_m_per_h[NDISC];

    int    n_sub_used;
} SoilFluxes;

typedef struct {
    double in_rain_m;         // infiltrated
    double excess_m;          // rejected
    double perc_m;
    double AET_m;
    double lateral_m;
    double delta_storage_m;   // \u03a3(\u03b8_out-\u03b8_in)*dz
    double residual_m;        // in - (outs) - \u0394S
} TimestepSoilVolumeBalance;

#ifdef __cplusplus
} // extern "C"
#endif

#endif // SOIL_TYPES_H
