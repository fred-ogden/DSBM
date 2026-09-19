#!/bin/tcsh
# Run the stateless Noah-MP soil-water formulation with full output

set force   = "10x_hourly_rainfall_PET.txt"
set outdir  = "output/noahmp_fc_limit"
set dz      = "0.1,0.3,0.6,1.0"
set zc      = "0.05,0.25,0.7,1.5"
set store0  = "0.34"         # target initial storage (m)
set verb    = "1"            # stdout verbosity
set tstart  = "0"            # starting index column
set timecol = "index"        # choose index | jd | ymdh

rm -rf "$outdir"
mkdir -p "$outdir"

bin/soil_driver \
  --solver noahmp \
  --apply-fc-perc-threshold \
  --forcing "$force" \
  --outdir "$outdir" \
  --verbosity $verb \
  --write-theta --write-fluxes --write-mb \
  --timecol $timecol --tindex-start $tstart \
  --dz "$dz" --zc "$zc" \
  --init-target-storage $store0

set run_status = $status

if ($run_status != 0) then
    echo "ERROR: Noah-MP run failed with status $run_status"
    exit $run_status
endif

echo "Noah-MP run complete.  Output written to $outdir"
