#!/bin/tcsh -f
#
# run_tapenade_equal_head.csh
#
# One-command build and run of the E11 experiment (FLO 11/06/2026): the
# Tapenade tangent of the equal-head limit on downward transfers
# (limit_downward_transfer_to_equal_head() in src/soil_helpers.c),
# compared with the analytic implicit-function derivative and with
# central finite differences at every limiter activation on the forcing
# record, for a list of parameter sets.  See equal_head_check_driver.c.
#
# Usage, from the repository root (or from tapenade/equal_head/ as ./run_tapenade_equal_head.csh):
#     tapenade/equal_head/run_tapenade_equal_head.csh
#
# Steps:
#   1. clean rebuild of the production library
#   2. stage the differentiated sources in tapenade_input/equal_head/
#   3. run Tapenade -tangent into generated_d/
#   4. print the generated prototype
#   5. compile the generated code, the adapter, the driver, and the DSBM
#      kernel with the limiter call renamed to the driver's recording
#      wrapper; link
#   6. run each parameter set; results go to output/tapenade_equal_head/
#
# Everything generated lives in tapenade_input/ or output/, both ignored
# by git.  ASCII only.

set nonomatch

# Work from the repository root, wherever this script is started from:
# this script lives in tapenade/equal_head/, 2 levels below the root.
set script_dir = $0:h
if ( "$script_dir" == "$0" ) set script_dir = .
cd $script_dir/../..
if ( ! -f Makefile || ! -d src || ! -d include ) then
    echo "ERROR: could not find the repository root from $0"
    exit 1
endif

# ---- tunable parameters and directories ----
set experiment_dir = tapenade/equal_head
set work_dir       = tapenade_input/equal_head
set gen_d          = $work_dir/generated_d
set results_dir    = output/tapenade_equal_head
set exe            = bin/tapenade_equal_head_check
set forcing        = forcing/rain_pet_example.csv
set ndisc          = 4
set head_spec      = "equal_head_check_primal(limited_transfer_m)/(transfer_requested_m theta_upper theta_lower theta_sat phi_sat_cm b_exp)"
set sources        = "equal_head_check_primal.c soil_helpers.c"
set cflags         = "-std=c11 -O2 -Wall -Wextra -DNDISC=$ndisc"
# parameter sets "k_lf perc_limiter Ksat b": the config values, then the
# E10-bound corners in which tools/run_equal_head_census.csh found the
# limiter on the Newton path (census sets 4, 5, 9, 13, 17)
set parameter_sets = ( \
    "0.000005 0.81 1.2168 4.05" \
    "1e-5 0.01 100 2" \
    "1e-5 0.01 100 15" \
    "1e-5 0.999 100 15" \
    "1 0.01 100 15" \
    "1 0.999 100 15" )

if ( ! -f $forcing ) then
    echo "ERROR: $forcing not found"
    exit 1
endif
which tapenade >& /dev/null
if ( $status != 0 ) then
    echo "ERROR: tapenade is not on PATH"
    exit 1
endif

echo "==== 1. building production library (clean rebuild)"
make -s veryclean
make -s
if ( $status != 0 ) then
    echo "ERROR: make failed"
    exit 1
endif

echo "==== 2. preparing Tapenade input in $work_dir"
rm -rf $work_dir
mkdir -p $gen_d
cp include/*.h $work_dir/
( echo "#define NDISC $ndisc" ; cat include/soil_config.h ) >! $work_dir/soil_config.h
foreach source_file ( $experiment_dir/equal_head_check_primal.c src/soil_helpers.c )
    ( echo "#define NDISC $ndisc" ; cat $source_file ) >! $work_dir/$source_file:t
end

echo "==== 3. running Tapenade (tangent mode)"
cd $work_dir
tapenade -tangent -head "$head_spec" -O generated_d $sources >& generated_d/tapenade.log
cd ../..
if ( ! -f $gen_d/equal_head_check_primal_d.c || ! -f $gen_d/soil_helpers_d.c ) then
    echo "ERROR: Tapenade did not produce the tangent sources.  Log:"
    cat $gen_d/tapenade.log
    exit 1
endif
echo "---- Tapenade log (last 15 lines):"
tail -15 $gen_d/tapenade.log

echo ""
echo "==== 4. generated prototype"
grep -A12 "void equal_head_check_primal_d *(" $gen_d/equal_head_check_primal_d.c | sed -n '1,/)/p'

echo ""
echo "==== 5. compiling"
cc $cflags -I$gen_d -I$work_dir -c $gen_d/soil_helpers_d.c -o $work_dir/soil_helpers_d.o
if ( $status != 0 ) then
    echo "ERROR: compiling $gen_d/soil_helpers_d.c failed"
    exit 1
endif
cc $cflags -I$gen_d -I$work_dir -I$experiment_dir -c $experiment_dir/equal_head_tangent_adapter.c -o $work_dir/tangent_adapter.o
if ( $status != 0 ) then
    echo "ERROR: compiling the tangent adapter failed (generated signature differs from expected)"
    exit 1
endif
cc $cflags -Iinclude -I$experiment_dir -c $experiment_dir/equal_head_check_driver.c -o $work_dir/equal_head_check_driver.o
if ( $status != 0 ) exit 1
# the DSBM kernel with its limiter call routed to the driver's recording wrapper
cc $cflags -Iinclude -Dlimit_downward_transfer_to_equal_head=equal_head_recording_limit -c src/dsbm_soilmoisture_stateless.c -o $work_dir/dsbm_kernel_recording.o
if ( $status != 0 ) exit 1
mkdir -p bin
cc -o $exe $work_dir/equal_head_check_driver.o $work_dir/dsbm_kernel_recording.o $work_dir/tangent_adapter.o $work_dir/soil_helpers_d.o build/lib/libsoil.a -lm
if ( $status != 0 ) then
    echo "ERROR: link failed"
    exit 1
endif

echo ""
echo "==== 6. running"
rm -rf $results_dir
mkdir -p $results_dir
set summary = $results_dir/summary.txt
echo "E11 equal-head limiter derivative check, full forcing record" >! $summary
set n = 0
set n_fail = 0
foreach parameter_set ( $parameter_sets:q )
    @ n = $n + 1
    set values = ( $parameter_set )
    set set_dir = $results_dir/set_$n
    mkdir -p $set_dir
    echo "---- set $n : k_lf $values[1]  perc_limiter $values[2]  Ksat $values[3]  b $values[4]"
    $exe $set_dir $forcing $values[1] $values[2] $values[3] $values[4] >! $set_dir/results.txt
    set run_status = $status
    if ( $run_status == 1 ) then
        echo "ERROR: $exe failed for set $n; see $set_dir/results.txt"
        exit 1
    endif
    if ( $run_status != 0 ) @ n_fail = $n_fail + 1
    cat $set_dir/results.txt
    echo "" >> $summary
    echo "set $n : $parameter_set" >> $summary
    grep -e "activations:" -e "Newton-path activations failing" -e "zero-transfer activations" -e "RESULT" $set_dir/results.txt >> $summary
end

echo ""
cat $summary
if ( $n_fail != 0 ) then
    echo "$n_fail parameter set(s) FAILED.  Summary: $summary"
else
    echo "All parameter sets PASSED.  Summary: $summary"
endif
echo "Done."
