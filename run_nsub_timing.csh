#!/bin/tcsh -f
#
# run_nsub_timing.csh
#
# Runs the full DSBM simulation (configs/soil_params.dat, all of
# forcing/rain_pet_example.csv) once for each substep setting and
# tabulates cost, wall-clock time, and the main volume-balance terms.
#
#   --nsub 0   original adaptive substep choice
#   --nsub N   fixed at N substeps per timestep (default is 4)
#
# Usage, from the repository root:
#     ./run_nsub_timing.csh                 (default list: 0 1 2 4 6 8 12 24)
#     ./run_nsub_timing.csh 0 4 12          (any list of settings)
#
# Each run writes to output/nsub_<N>/.  The table is also saved to
# output/nsub_timing.txt.  Wall-clock time is as reported by soil_driver
# (simulation loop only); repeat the script to judge timing noise.
#
# ASCII only.

set nonomatch
set settings = ( 0 1 2 4 6 8 12 24 )
if ( $#argv > 0 ) set settings = ( $argv )

set force  = "forcing/rain_pet_example.csv"
set dz     = "0.1,0.3,0.6,1.0"
set zc     = "0.05,0.25,0.7,1.5"
set table  = output/nsub_timing.txt

make -s
if ( $status != 0 ) then
    echo "ERROR: make failed"
    exit 1
endif
mkdir -p output

# header
printf "%-8s %12s %9s %10s %12s %12s %12s %12s\n" \
    "nsub" "substeps" "mean/h" "wall_s" "lateral_m" "perc_m" "AET_m" "residual_m" >! $table
printf "%-8s %12s %9s %10s %12s %12s %12s %12s\n" \
    "--------" "------------" "---------" "----------" "------------" "------------" "------------" "------------" >> $table

foreach n ( $settings )
    set outdir = output/nsub_$n
    rm -rf $outdir
    mkdir -p $outdir

    bin/soil_driver \
        --config configs/soil_params.dat \
        --solver dsbm \
        --forcing $force \
        --outdir $outdir \
        --verbosity 0 \
        --write-volbal \
        --timestamp datetime \
        --dz "$dz" --zc "$zc" \
        --nsub $n
    if ( $status != 0 ) then
        echo "ERROR: soil_driver failed for --nsub $n"
        exit 1
    endif

    set summary  = $outdir/volbal_summary.out
    set substeps = `grep "Total number of sub-timesteps" $summary | awk '{print $NF}'`
    set mean_h   = `grep "Mean sub-timesteps per output step" $summary | awk '{print $NF}'`
    set wall_s   = `grep "wall-clock" $summary | awk '{print $(NF-1)}'`
    set lateral  = `grep "Lateral flow generated" $summary | awk '{print $(NF-1)}'`
    set perc     = `grep "Percolation to g.w." $summary | awk '{print $(NF-1)}'`
    set aet      = `grep "Evapotranspiration" $summary | awk '{print $(NF-1)}'`
    set resid    = `grep "Volume balance residual  " $summary | awk '{print $(NF-1)}'`

    set label = "$n"
    if ( "$n" == "0" ) set label = "0 adapt"

    printf "%-8s %12s %9s %10s %12s %12s %12s %12s\n" \
        "$label" "$substeps" "$mean_h" "$wall_s" "$lateral" "$perc" "$aet" "$resid" >> $table
end

cat $table
echo ""
echo "Table saved to $table"
