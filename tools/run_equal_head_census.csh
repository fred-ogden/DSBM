#!/bin/tcsh -f
#
# run_equal_head_census.csh
#
# How often the equal-head limit on downward transfers is active, over
# the full forcing record, for a list of soil parameter sets
# (FLO 11/06/2026).  Uses the census build (bin/soil_driver_census) and
# reports, per parameter set, the census lines for
#   eq_head_limited      downward transfer reduced by the limit
#   eq_head_zero         limit set the downward transfer to 0
#   eq_head_late_bisect  bisection step in the last 5 limiter iterations
# The last one matters for differentiability: Tapenade's derivative of
# the limited transfer is the correct implicit-function derivative only
# if the final iterations are Newton steps.
#
# Each parameter set is "k_lf perc_limiter Ksat b" (m/h, -, cm/h, -).
# Each run uses a copy of configs/soil_params.dat with those four values
# replaced; phi_sat follows Ksat through the config's
# soil_saturated_capillary_head_calc_from_ksat=TRUE, and theta_fc and
# theta_wp follow phi_sat and b.  Everything else is as in the config.
#
# The default list is the config values plus the 16 corners of the
# bounds clamped by the E10 BFGS search
# (tapenade/soil_params/tapenade_soil_params_driver.c):
#   k_lf 1e-5 to 1 m/h, perc_limiter 0.01 to 0.999,
#   Ksat 0.01 to 100 cm/h, b 2 to 15.
# Replace the list with your calibration bounds as needed.
#
# Usage, from the repository root (or from tools/ as ./run_equal_head_census.csh):
#     tools/run_equal_head_census.csh
#
# Output: output/equal_head_census/census_<n>.txt for each parameter
# set, and output/equal_head_census/summary.txt.  ASCII only.

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

# ---- tunable parameters and directories ----
set results_dir = output/equal_head_census
set base_config = configs/soil_params.dat
set force       = "forcing/rain_pet_example.csv"
set dz          = "0.1,0.3,0.6,1.0"
set zc          = "0.05,0.25,0.7,1.5"
# parameter sets: "k_lf perc_limiter Ksat b"
set parameter_sets = ( \
    "0.000005 0.81 1.2168 4.05" \
    "1e-5 0.01 0.01 2"    "1e-5 0.01 0.01 15"    "1e-5 0.01 100 2"    "1e-5 0.01 100 15" \
    "1e-5 0.999 0.01 2"   "1e-5 0.999 0.01 15"   "1e-5 0.999 100 2"   "1e-5 0.999 100 15" \
    "1 0.01 0.01 2"       "1 0.01 0.01 15"       "1 0.01 100 2"       "1 0.01 100 15" \
    "1 0.999 0.01 2"      "1 0.999 0.01 15"      "1 0.999 100 2"      "1 0.999 100 15" )

set klf_key  = "soil_reservoir_rate_const_to_subsurface_lateral_flow"
set perc_key = "soil_to_gw_percolation_rate_limiter_0_to_1"
set ksat_key = "soil_sat_hydraulic_conductivity_cm_per_h"
set b_key    = "soil_Clapp_Hornberger_exponent_b"

if ( ! -f $base_config ) then
    echo "ERROR: $base_config not found"
    exit 1
endif
if ( ! -f $force ) then
    echo "ERROR: $force not found"
    exit 1
endif

make -s census
if ( $status != 0 ) then
    echo "ERROR: make census failed"
    exit 1
endif
rm -rf $results_dir
mkdir -p $results_dir

set summary = $results_dir/summary.txt
echo "Equal-head limit census, full forcing record" >! $summary
echo "set  k_lf(m/h)  perc_limiter  Ksat(cm/h)  b" >> $summary

set n = 0
foreach parameter_set ( $parameter_sets:q )
    @ n = $n + 1
    set values = ( $parameter_set )
    set config = $results_dir/soil_params_$n.dat
    sed -e "s/^${klf_key}=[^[]*/${klf_key}=$values[1]/" -e "s/^${perc_key}=[^[]*/${perc_key}=$values[2]/" -e "s/^${ksat_key}=[^[]*/${ksat_key}=$values[3]/" -e "s/^${b_key}=[^[]*/${b_key}=$values[4]/" $base_config >! $config

    set outdir = $results_dir/run_$n
    set report = $results_dir/census_$n.txt
    mkdir -p $outdir

    echo "==== set $n : k_lf $values[1]  perc_limiter $values[2]  Ksat $values[3]  b $values[4]"
    bin/soil_driver_census --config $config --solver dsbm --forcing $force --outdir $outdir --verbosity 0 --write-volbal --timestamp datetime --dz "$dz" --zc "$zc" >! $report
    if ( $status != 0 ) then
        echo "ERROR: soil_driver_census failed for set $n"
        exit 1
    endif

    echo "" >> $summary
    echo "set $n : $parameter_set" >> $summary
    grep "residual  " $outdir/volbal_summary.out >> $summary
    grep -A 3 "^eq_head_limited" $report >> $summary
    grep -A 3 "^eq_head_zero" $report >> $summary
    grep -A 1 "^eq_head_late_bisect" $report >> $summary
    grep -A 1 "^eq_head_call" $report >> $summary
    grep -A 1 "^eq_head_prescreen" $report >> $summary
end

cat $summary
echo "Summary written to $summary"
echo "Done."
