#!/bin/tcsh
# Run the stateless DSBM soil column driver with full output

set force   = "forcing/rain_pet_example.csv"
set outdir  = "output/dsbm"
set dz      = "0.1,0.3,0.6,1.0"
set zc      = "0.05,0.25,0.7,1.5"
set verb    = "1"            # stdout verbosity
set timecol = "datetime"     # choose timestep | datetime | juliandate

rm -rf "$outdir"
mkdir -p "$outdir"

bin/soil_driver \
  --config configs/soil_params.dat \
  --solver dsbm \
  --forcing "$force" \
  --outdir "$outdir" \
  --verbosity $verb \
  --write-theta --write-fluxes --write-volbal \
  --timestamp $timecol \
  --use-lut --lut-n 400 --lut-Theta-min 1e-6 \
  --dz "$dz" --zc "$zc" \


set run_status = $status

if ($run_status != 0) then
    echo "ERROR: DSBM run failed with status $run_status"
    exit $run_status
endif

echo "DSBM run complete.  Output written to $outdir"
