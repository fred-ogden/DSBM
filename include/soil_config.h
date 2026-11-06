/*
 * soil_config.h
 *
 * Compile-time configuration constants for the soil-moisture calculation modules.
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

#ifndef SOIL_CONFIG_H
#define SOIL_CONFIG_H

// Compile-time discretization count.
#ifndef NDISC
#error "NDISC must be defined at compile time (number of soil discs)."
#endif

// The one dry-end limit on soil moisture (FLO 11/06/2026).
//
// The Clapp-Hornberger (CH) relations are written in the scaled soil
// moisture THETA = theta / theta_e, where theta_e is the effective
// porosity (theta_sat in the code).  They assume theta_r = 0, so psi and
// K are singular at theta = 0 and need a lower limit.  From a soil
// physics standpoint THETA = 0.03 is a reasonable one, and it is the
// single limit used everywhere soil moisture is bounded below:
//   - the clamps in psi_from_theta(), K_from_theta() and theta_from_psi(),
//   - the minimum of the CH lookup table (which has 2 points, since the
//     CH relations are log-linear in THETA),
//   - the kernel's numerical floor theta_floor = 0.03 theta_e (or theta_r
//     when theta_r > 0: a retention relation with a residual soil
//     moisture is not CH, and its floor is theta_r),
//   - the derived wilting point, clamped at no less than 0.03 theta_e so
//     that AET (which stops at theta_wp) cannot take a disc below the
//     floor.
// It is not expected to be reached: AET stops at the wilting point.
// This replaces the absolute THETA_MIN = 1.0e-3 m3/m3 used before, which
// differed from the lookup-table minimum.
#define CH_SCALED_SOIL_MOISTURE_MIN 0.03   // (-), minimum of theta / theta_e

// Tiny epsilon to avoid divide-by-zero.
#ifndef SOIL_EPS
#define SOIL_EPS 1.0e-12
#endif

// Inline helper macro
#ifndef SOIL_INLINE
#define SOIL_INLINE static inline
#endif

#endif // SOIL_CONFIG_H
