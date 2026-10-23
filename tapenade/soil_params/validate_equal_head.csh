#!/bin/tcsh -f
# FLO 10/26: Rebuild and validate full-model and active-limiter derivatives.
set nonomatch
set script_dir = $0:h
if ( "$script_dir" == "$0" ) set script_dir = .
cd $script_dir/../..
if ( ! -f Makefile || ! -d tapenade/soil_params ) then
    echo "FAIL: cannot locate DSBM repository root"
    exit 1
endif
set experiment_dir = tapenade/soil_params
set work_dir = tapenade_input/equal_head
set results_dir = output/equal_head_validation
set kit_dir = /user2/ogden/opt/tapenade/tapenade_3.16/ADFirstAidKit
set flags = '-std=c99 -O2 -Wall -Wextra -DNDISC=4 -DTHETA_MIN=1.0e-3'
if ( ! -f $kit_dir/adStack.c ) then
    echo "FAIL: Tapenade ADFirstAidKit not found at $kit_dir"
    exit 1
endif
mkdir -p $results_dir

echo '==== Full DSBM cost-function differentiation (E10a, E10b)'
$experiment_dir/run_tapenade_soil_params.csh ab >& $results_dir/full_model_build.log
if ( $status != 0 ) then
    echo "FAIL: full-model build or execution; see $results_dir/full_model_build.log"
    exit 1
endif
python3 $experiment_dir/check_equal_head_validation.py full output/tapenade_soil_params/results.txt
if ( $status != 0 ) exit 1

echo '==== Active equal-head limiter: generating tangent and adjoint'
rm -rf $work_dir
mkdir -p $work_dir/generated_d $work_dir/generated_b
cp include/*.h $work_dir/
( echo '#define NDISC 4' ; cat $experiment_dir/equal_head_ad_primal.c ) >! $work_dir/equal_head_ad_primal.c
( echo '#define NDISC 4' ; cat src/soil_helpers.c ) >! $work_dir/soil_helpers.c
cd $work_dir
tapenade -tangent -head 'equal_head_ad_primal(limited_transfer_m)/(b_exp)' -O generated_d equal_head_ad_primal.c soil_helpers.c >& generated_d/tapenade.log
if ( $status != 0 || ! -f generated_d/equal_head_ad_primal_d.c ) then
    echo 'FAIL: tangent generation; inspect generated_d/tapenade.log'
    exit 1
endif
tapenade -reverse -head 'equal_head_ad_primal(limited_transfer_m)/(b_exp)' -O generated_b equal_head_ad_primal.c soil_helpers.c >& generated_b/tapenade.log
if ( $status != 0 || ! -f generated_b/equal_head_ad_primal_b.c ) then
    echo 'FAIL: adjoint generation; inspect generated_b/tapenade.log'
    exit 1
endif
cd ../..

echo '==== Compiling active-limiter test drivers'
cc $flags -I$work_dir -I$work_dir/generated_d $experiment_dir/test_equal_head_derivative.c $work_dir/generated_d/equal_head_ad_primal_d.c $work_dir/generated_d/soil_helpers_d.c build/lib/libsoil.a -lm -o $work_dir/test_equal_head_derivative
if ( $status != 0 ) exit 1
cc $flags -D_POSIX_C_SOURCE=200809L -I$work_dir -I$work_dir/generated_b -I$kit_dir $experiment_dir/test_equal_head_adjoint.c $work_dir/generated_b/equal_head_ad_primal_b.c $work_dir/generated_b/soil_helpers_b.c $kit_dir/adStack.c build/lib/libsoil.a -lm -o $work_dir/test_equal_head_adjoint
if ( $status != 0 ) exit 1

echo '==== Running active-limiter tests'
$work_dir/test_equal_head_derivative >! $results_dir/tangent.txt
if ( $status != 0 ) exit 1
$work_dir/test_equal_head_adjoint >! $results_dir/adjoint.txt
if ( $status != 0 ) exit 1
python3 $experiment_dir/check_equal_head_validation.py active $results_dir/tangent.txt $results_dir/adjoint.txt
if ( $status != 0 ) exit 1
echo 'PASS: all DSBM equal-head differentiation checks'
echo "Results: $results_dir"
