#!/bin/tcsh
# Run the stateless Noah-MP soil-water formulation with full output

set force   = "forcing/rain_pet_example.csv"
set outdir  = "output/noahmp_fc_limit"
set dz      = "0.1,0.3,0.6,1.0"
set zc      = "0.05,0.25,0.7,1.5"
set verb    = "1"            # stdout verbosity
set timecol = "datetime"     # choose timestep | datetime | juliandate

rm -rf "$outdir"
mkdir -p "$outdir"

bin/soil_driver \
  --config configs/soil_params.dat \
  --solver noahmp \
  --apply-fc-perc-threshold \
  --forcing "$force" \
  --outdir "$outdir" \
  --verbosity $verb \
  --write-theta --write-fluxes --write-volbal \
  --timestamp $timecol \
  --dz "$dz" --zc "$zc" \

set run_status = $status

if ($run_status != 0) then
    echo "ERROR: Noah-MP run failed with status $run_status"
    exit $run_status
endif

echo "Noah-MP run complete.  Output written to $outdir"
