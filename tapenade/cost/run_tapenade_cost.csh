#!/bin/tcsh -f
#
# run_tapenade_cost.csh
#
# One-command build and run of the E8 cost-function gradient experiments:
# Tapenade tangent of a scalar cost function (sum of 1 - NSE for hourly
# lateral flow and percolation against synthetic observations) with
# respect to k_lf and perc_limiter.
#   E8a gradient check, E8b one-parameter sweeps, E8c BFGS calibration.
#
# Usage, from the repository root:
#     ./run_tapenade_cost.csh           (E8a, E8b, E8c)
#     ./run_tapenade_cost.csh c         (E8c only; any letter string works)
#
# Steps:
#   1. build the production library (make)
#   2. copy the three differentiated sources and the headers into
#      tapenade_input/cost/, with NDISC defined in the copied soil_config.h
#      (so Tapenade needs no -D option)
#   3. run Tapenade in tangent mode on dsbm_cost_from_params
#   4. print the generated tangent prototype
#   5. compile the generated code, the adapter, the primal wrapper and
#      the test driver, and link
#   6. run the tests; results go to output/tapenade_cost/
#
# Everything generated lives in tapenade_input/ or output/, both ignored
# by git.  ASCII only.

set nonomatch
set experiments = abc
if ( $#argv > 0 ) set experiments = "$1"
set ndisc = 4
set theta_min = 1.0e-03
set work_dir = tapenade_input/cost
set gen_dir = $work_dir/generated
set results_dir = output/tapenade_cost
set exe = tapenade_cost_derivative_test
set head_function = dsbm_cost_from_params
set dependents = "cost_function_value"
set independents = "klf_m_per_h perc_limiter_0_to_1"

set cflags = "-std=c11 -O2 -Wall -Wextra -DNDISC=$ndisc -DTHETA_MIN=$theta_min"

# Tapenade support headers, if present (harmless if not needed)
set tapenade_home = /user2/ogden/opt/tapenade/tapenade_3.16
set tapenade_kit_flags = ""
if ( -d $tapenade_home/ADFirstAidKit ) then
    set tapenade_kit_flags = "-I$tapenade_home/ADFirstAidKit"
endif

echo "==== 1. building production library (clean rebuild)"
# Always rebuild from scratch: a stale libsoil.a compiled against an older
# SoilControl layout silently corrupts every primal (finite-difference) run.
make -s veryclean
make -s
if ( $status != 0 ) then
    echo "ERROR: make failed"
    exit 1
endif

echo "==== 2. preparing Tapenade input in $work_dir"
rm -rf $work_dir
mkdir -p $gen_dir
cp include/*.h $work_dir/
( echo "#define NDISC $ndisc" ; cat include/soil_config.h ) >! $work_dir/soil_config.h
# soil_helpers.c tests NDISC on its first lines, before any #include,
# so the definition must be at the top of each copied source file too.
foreach source_file ( tapenade_cost_primal.c src/dsbm_soilmoisture_stateless.c src/soil_helpers.c )
    ( echo "#define NDISC $ndisc" ; cat $source_file ) >! $work_dir/$source_file:t
end

echo "==== 3. running Tapenade (tangent mode)"
cd $work_dir
tapenade -tangent \
    -head "${head_function}(${dependents})/(${independents})" \
    -O generated \
    tapenade_cost_primal.c dsbm_soilmoisture_stateless.c soil_helpers.c \
    >& generated/tapenade.log
cp generated/tapenade.log generated/tapenade_head_form.log

# find the generated file that defines the head's tangent
set head_file = `grep -ls "${head_function}_d *(" generated/*_d.c`
if ( "$head_file" == "" ) then
    echo "NOTE: -head 'f(outs)/(ins)' form produced no tangent; retrying with -vars/-outvars"
    rm -f generated/*_d.c generated/*_d.h generated/*.msg
    tapenade -tangent \
        -head ${head_function} \
        -vars "${independents}" \
        -outvars "${dependents}" \
        -O generated \
        tapenade_cost_primal.c dsbm_soilmoisture_stateless.c soil_helpers.c \
        >& generated/tapenade.log
    set head_file = `grep -ls "${head_function}_d *(" generated/*_d.c`
endif
cd ../..

if ( "$head_file" == "" ) then
    echo "ERROR: Tapenade did not produce ${head_function}_d."
    echo "---- log of first attempt (-head form):"
    cat $gen_dir/tapenade_head_form.log
    echo "---- log of second attempt (-vars/-outvars form):"
    cat $gen_dir/tapenade.log
    exit 1
endif

echo "Tapenade log (last 25 lines):"
tail -25 $gen_dir/tapenade.log
echo ""
echo "Generated files:"
ls -1 $gen_dir

# The adapter #includes tapenade_cost_primal_d.c.  If Tapenade put the
# head's tangent in a differently named file, copy it to that name.
if ( "$head_file" != "generated/tapenade_cost_primal_d.c" ) then
    echo "NOTE: head tangent is in $head_file; copying to tapenade_cost_primal_d.c"
    cp $work_dir/$head_file $gen_dir/tapenade_cost_primal_d.c
endif

echo ""
echo "==== 4. generated tangent prototype (compare with tapenade_cost_tangent_adapter.c)"
grep -A12 "void ${head_function}_d *(" $gen_dir/tapenade_cost_primal_d.c | sed -n '1,/)/p'

# Does the generated prototype carry an integer tangent slot for
# n_sub_setting?  The adapter adapts at compile time.
set nsub_slot_flag = ""
set nsub_slot_count = `grep -A12 "void ${head_function}_d *(" $gen_dir/tapenade_cost_primal_d.c | sed -n '1,/)/p' | grep -c n_sub_settingd`
if ( "$nsub_slot_count" != "0" ) then
    set nsub_slot_flag = "-DCOST_TANGENT_HAS_NSUB_SLOT"
    echo "(prototype has an n_sub_settingd slot; adapter passes 0 for it)"
else
    echo "(prototype has no n_sub_settingd slot)"
endif

echo ""
echo "==== 5. compiling"
set objects = ""
foreach gen_c ( $gen_dir/*_d.c )
    # the head file is compiled through the adapter, not separately
    if ( "$gen_c" == "$gen_dir/tapenade_cost_primal_d.c" ) continue
    if ( "$gen_c" == "$work_dir/$head_file" ) continue
    set obj = $gen_c:r.o
    cc $cflags -I$gen_dir -I$work_dir $tapenade_kit_flags -c $gen_c -o $obj
    if ( $status != 0 ) then
        echo "ERROR: compiling $gen_c failed"
        exit 1
    endif
    set objects = "$objects $obj"
end

cc $cflags $nsub_slot_flag -I$gen_dir -I$work_dir $tapenade_kit_flags \
    -c tapenade_cost_tangent_adapter.c -o $gen_dir/tapenade_cost_tangent_adapter.o
if ( $status != 0 ) then
    echo "ERROR: compiling the tangent adapter failed (generated signature differs from expected)"
    exit 1
endif

cc $cflags -Iinclude -c tapenade_cost_primal.c -o $gen_dir/tapenade_cost_primal.o
if ( $status != 0 ) exit 1
cc $cflags -Iinclude -c tapenade_cost_driver.c -o $gen_dir/tapenade_cost_driver.o
if ( $status != 0 ) exit 1

set link_objects = "$gen_dir/tapenade_cost_driver.o $gen_dir/tapenade_cost_primal.o $gen_dir/tapenade_cost_tangent_adapter.o $objects"
cc -o $exe $link_objects build/lib/libsoil.a -lm
if ( $status != 0 ) then
    # Tapenade sometimes emits copies of primal routines; they are the
    # same code as libsoil.a, so allow the duplicate definitions.
    echo "NOTE: retrying link with --allow-multiple-definition"
    cc -o $exe $link_objects build/lib/libsoil.a -lm -Wl,--allow-multiple-definition
    if ( $status != 0 ) then
        echo "ERROR: link failed"
        exit 1
    endif
endif

echo ""
echo "==== 6. running experiments"
rm -rf $results_dir
mkdir -p $results_dir
./$exe $results_dir forcing/rain_pet_example.csv $experiments | tee $results_dir/results.txt
set run_status = $status

echo ""
echo "Results: $results_dir/results.txt and $results_dir/*.csv"
exit $run_status
