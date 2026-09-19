#ifndef SOIL_KERNEL_STATELESS_H
#define SOIL_KERNEL_STATELESS_H

#include "soil_data_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int soil_step_one_hour_stateless(
    const SoilControl        *ctrl,
    const SoilGeometry       *geom,
    const SoilParameters     *par,
    const SoilLookupTables   *lut,      // may be NULL to use analytic CH
    const SoilStateIn        *sin,
    const SoilForcing        *forcing,
    SoilStateOut             *sout,
    SoilFluxes               *flux,
    TimestepSoilMassbal      *mb,
    FILE                     *debug_fptr);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // SOIL_KERNEL_STATELESS_H
