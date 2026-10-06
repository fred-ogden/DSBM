/*
 * noahmp_soilmoisture_stateless.h
 *
 * Public interface for the stateless Noah-MP comparison soil-moisture solver.
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

#ifndef NOAHMP_SOILWATER_STATELESS_H
#define NOAHMP_SOILWATER_STATELESS_H

#include <stdio.h>
#include "soil_data_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int noahmp_soil_step_one_hour_stateless(
    const SoilControl        *ctrl,
    const SoilGeometry       *geom,
    const SoilParameters     *par,
    const SoilStateIn        *sin,
    const SoilForcing        *forcing,
    SoilStateOut             *sout,
    SoilFluxes               *flux,
    TimestepSoilVolumeBalance      *volbal,
    FILE                     *debug_fptr);

#ifdef __cplusplus
}
#endif

#endif
