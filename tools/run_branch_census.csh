#!/bin/tcsh -f
#
# run_branch_census.csh
#
# Branch-activity census of the DSBM kernel over the full forcing record
# (configs/soil_params.dat, forcing/rain_pet_example.csv).  Counts how
# often each threshold, cap and branch is active: Darcy-Buckingham flux
# direction and transfer caps, acceptance-chain caps, percolation cap,
# ET stress regimes and cap, lateral activity, theta_fc and theta_wp
# crossings, and roundoff clamps.
#
# Runs once per lateral rate constant k_lf (m/h).  Each run uses a copy
# of the config with only k_lf replaced; all other parameters are as in
# configs/soil_params.dat (including perc_limiter = 0.81).
#
# Usage, from the repository root (or from tools/ as ./run_branch_census.csh):
#     tools/run_branch_census.csh                      (k_lf = config, 1e-3, 1e-2)
#     tools/run_branch_census.csh 1e-3                 (any list of k_lf values)
#     tools/run_branch_census.csh --nsub 0 1e-3        (extra soil_driver options first)
#
# Output: output/census/census_klf_<value>.txt for each k_lf.
#
# ASCII only.

set nonomatch

# Work from the repository root, wherever this script is started from:
# this script lives in tools/, 1 level below the root.
set script_dir = $0:h
if ( "$script_dir" == "$0" ) set script_dir = .
cd $script_dir/..
if ( ! -f Makefile || ! -d src || ! -d include ) then
    echo "ERROR: could not find the repository root from $0"
    exit 1
endif
set driver_options = ( )
set klf_values = ( )

# options start with "--" and take one value except flag-only options
while ( $#argv > 0 )
    if ( "$1" == "--nsub" || "$1" == "--n-sub-minimum" ) then
        set driver_options = ( $driver_options $1 $2 )
        shift ; shift
    else if ( "$1" =~ --* ) then
        set driver_options = ( $driver_options $1 )
        shift
    else
        set klf_values = ( $klf_values $1 )
        shift
    endif
end
if ( $#klf_values == 0 ) set klf_values = ( config 1e-3 1e-2 )

set force  = "forcing/rain_pet_example.csv"
set dz     = "0.1,0.3,0.6,1.0"
set zc     = "0.05,0.25,0.7,1.5"
set klf_key = "soil_reservoir_rate_const_to_subsurface_lateral_flow"

make -s census
if ( $status != 0 ) then
    echo "ERROR: make census failed"
    exit 1
endif
mkdir -p output/census

foreach klf ( $klf_values )
    set config = configs/soil_params.dat
    if ( "$klf" != "config" ) then
        set config = output/census/soil_params_klf_$klf.dat
        sed "s/^${klf_key}=[^[]*/${klf_key}=${klf}/" configs/soil_params.dat >! $config
    endif

    set outdir = output/census/run_klf_$klf
    set report = output/census/census_klf_$klf.txt
    rm -rf $outdir
    mkdir -p $outdir

    echo "==== k_lf = $klf   options: $driver_options"
    grep "^${klf_key}" $config

    bin/soil_driver_census \
        --config $config \
        --solver dsbm \
        --forcing $force \
        --outdir $outdir \
        --verbosity 0 \
        --write-volbal \
        --timestamp datetime \
        --dz "$dz" --zc "$zc" \
        $driver_options >! $report
    if ( $status != 0 ) then
        echo "ERROR: soil_driver_census failed for k_lf = $klf"
        exit 1
    endif
    grep "residual  " $outdir/volbal_summary.out
    echo "census written to $report"
end
