#!/bin/tcsh
# Run the stateless DSBM soil column driver with full output

set force   = "10x_hourly_rainfall_PET.txt"
set outdir  = "output/dsbm_fc_limit"
set dz      = "0.1,0.3,0.6,1.0"
set zc      = "0.05,0.25,0.7,1.5"
set store0  = "0.34"         # target initial storage (m)
set verb    = "1"            # stdout verbosity
set tstart  = "0"            # starting index column
set timecol = "index"        # choose index | jd | ymdh

rm -rf "$outdir"
mkdir -p "$outdir"

bin/soil_driver \
  --solver dsbm \
  --apply-fc-perc-threshold \
  --forcing "$force" \
  --outdir "$outdir" \
  --verbosity $verb \
  --write-theta --write-fluxes --write-mb \
  --timecol $timecol --tindex-start $tstart \
  --use-lut --lut-n 400 --lut-Theta-min 1e-6 \
  --dz "$dz" --zc "$zc" \
  --init-target-storage $store0

set run_status = $status

if ($run_status != 0) then
    echo "ERROR: DSBM run failed with status $run_status"
    exit $run_status
endif

echo "DSBM run complete.  Output written to $outdir"
