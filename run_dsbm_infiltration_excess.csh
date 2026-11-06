#!/bin/tcsh
# Run the stateless DSBM soil driver with full output, with
# control_soil_infiltration_flux_upper_boundary_true_false=FALSE (rain enters
# disc 1 only; the rest is infiltration excess).  DSBM's original form, for
# testing; TRUE (run_dsbm.csh) is required for equivalence to WRF-Hydro Noah-MP.
# FLO 11/06/2026.

set force   = "forcing/rain_pet_example.csv"
set outdir  = "output/dsbm_infiltration_excess"
set dz      = "0.1,0.3,0.6,1.0"
set zc      = "0.05,0.25,0.7,1.5"
set verb    = "1"            # stdout verbosity
set timecol = "datetime"     # choose timestep | datetime | juliandate

rm -rf "$outdir"
mkdir -p "$outdir"

bin/soil_driver \
  --config configs/soil_params_infiltration_excess.dat \
  --solver dsbm \
  --forcing "$force" \
  --outdir "$outdir" \
  --verbosity $verb \
  --write-theta --write-fluxes --write-volbal \
  --timestamp $timecol \
  --use-lut \
  --dz "$dz" --zc "$zc" \


set run_status = $status

if ($run_status != 0) then
    echo "ERROR: DSBM infiltration-excess run failed with status $run_status"
    exit $run_status
endif

echo "DSBM infiltration-excess run complete.  Output written to $outdir"
