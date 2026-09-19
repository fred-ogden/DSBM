#ifndef NOAHMP_SOILWATER_STATELESS_H
#define NOAHMP_SOILWATER_STATELESS_H

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
    TimestepSoilMassbal      *mb,
    FILE                     *debug_fptr);

#ifdef __cplusplus
}
#endif

#endif
