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

// Fictitious lower bound on water content (m3/m3).
#ifndef THETA_MIN
#define THETA_MIN 1.0e-03
#endif

// Tiny epsilon to avoid divide-by-zero.
#ifndef SOIL_EPS
#define SOIL_EPS 1.0e-12
#endif

// Inline helper macro
#ifndef SOIL_INLINE
#define SOIL_INLINE static inline
#endif

#endif // SOIL_CONFIG_H
