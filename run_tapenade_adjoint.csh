#!/bin/tcsh -f
#
# run_tapenade_adjoint.csh
#
# One-command build and run of the E9 adjoint experiments: Tapenade
# REVERSE mode (adjoint) on the E8 cost-function wrapper
# dsbm_cost_from_params(), compared with the TANGENT mode from the same
# source.
#   E9a adjoint versus tangent and dot-product test, E9b run time,
#   E9c adjoint memory, E9d BFGS calibration with each gradient.
#
# Usage, from the repository root:
#     ./run_tapenade_adjoint.csh            (E9a-E9d)
#     ./run_tapenade_adjoint.csh ab         (any letter string)
#
# Steps:
#   1. clean rebuild of the production library
#   2. stage the differentiated sources in tapenade_input/adjoint/
#   3. run Tapenade twice: -tangent into generated_d/, -reverse into
#      generated_b/
#   4. print both generated prototypes; detect integer derivative slots
#   5. compile generated code, Tapenade's ADFirstAidKit stack (adStack.c,
#      needed by reverse mode), both adapters, the primal wrapper and the
#      driver; link
#   6. run; results go to output/tapenade_adjoint/
#
# The tangent and reverse files each carry identical copies of some
# undifferentiated helper routines, so the link allows duplicate
# definitions.  Everything generated lives in tapenade_input/ or output/,
# both ignored by git.  ASCII only.

set nonomatch
set experiments = abcd
if ( $#argv > 0 ) set experiments = "$1"
set ndisc = 4
set theta_min = 1.0e-03
set work_dir = tapenade_input/adjoint
set gen_d = $work_dir/generated_d
set gen_b = $work_dir/generated_b
set results_dir = output/tapenade_adjoint
set exe = tapenade_adjoint_test
set head_function = dsbm_cost_from_params
set head_spec = "dsbm_cost_from_params(cost_function_value)/(klf_m_per_h perc_limiter_0_to_1)"
set sources = "tapenade_cost_primal.c dsbm_soilmoisture_stateless.c soil_helpers.c"

set cflags = "-std=c11 -O2 -Wall -Wextra -DNDISC=$ndisc -DTHETA_MIN=$theta_min"

set tapenade_home = /user2/ogden/opt/tapenade/tapenade_3.16
set kit_dir = $tapenade_home/ADFirstAidKit
if ( ! -f $kit_dir/adStack.c || ! -f $kit_dir/adStack.h ) then
    echo "ERROR: reverse mode needs $kit_dir/adStack.c and adStack.h"
    echo "       (Tapenade's ADFirstAidKit).  Contents of $tapenade_home :"
    ls $tapenade_home
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
mkdir -p $gen_d $gen_b
cp include/*.h $work_dir/
( echo "#define NDISC $ndisc" ; cat include/soil_config.h ) >! $work_dir/soil_config.h
foreach source_file ( tapenade_cost_primal.c src/dsbm_soilmoisture_stateless.c src/soil_helpers.c )
    ( echo "#define NDISC $ndisc" ; cat $source_file ) >! $work_dir/$source_file:t
end

echo "==== 3. running Tapenade (tangent mode, then reverse mode)"
cd $work_dir
tapenade -tangent -head "$head_spec" -O generated_d $sources >& generated_d/tapenade.log
# Reverse mode.  Tapenade 3.16 crashed in its TBR recomputation analysis
# (ADTBRAnalyzer RecompInfo, ArrayIndexOutOfBounds) on this code, so retry
# with that optimization off, and finally with TBR analysis off entirely
# (stores every overwritten value: more memory, simplest for Tapenade).
set reverse_options_used = "none"
foreach reverse_options ( "" "-nooptim recomputeintermediates" "-nooptim recomputeintermediates -nooptim tbr" )
    rm -f generated_b/*
    tapenade -reverse $reverse_options -head "$head_spec" -O generated_b $sources >& generated_b/tapenade.log
    if ( -f generated_b/tapenade_cost_primal_b.c ) then
        set reverse_options_used = "$reverse_options"
        if ( "$reverse_options_used" == "" ) set reverse_options_used = "(defaults)"
        break
    endif
    echo "NOTE: tapenade -reverse $reverse_options failed; log saved"
    cp generated_b/tapenade.log failed_reverse_`echo "$reverse_options" | tr -c 'a-z' '_'`.log
end
cd ../..
echo "reverse-mode options that worked: $reverse_options_used"

foreach mode ( d b )
    set gen = $work_dir/generated_$mode
    if ( ! -f $gen/tapenade_cost_primal_$mode.c ) then
        echo "ERROR: Tapenade did not produce $gen/tapenade_cost_primal_$mode.c.  Log:"
        cat $gen/tapenade.log
        echo "---- logs of earlier reverse attempts:"
        cat $work_dir/failed_reverse_*.log
        exit 1
    endif
    echo "---- Tapenade log, mode $mode (last 15 lines):"
    tail -15 $gen/tapenade.log
end

echo ""
echo "==== 4. generated prototypes"
grep -A12 "void ${head_function}_d *(" $gen_d/tapenade_cost_primal_d.c | sed -n '1,/)/p'
grep -A12 "void ${head_function}_b *(" $gen_b/tapenade_cost_primal_b.c | sed -n '1,/)/p'

set tangent_flag = ""
set count = `grep -A12 "void ${head_function}_d *(" $gen_d/tapenade_cost_primal_d.c | sed -n '1,/)/p' | grep -c n_sub_settingd`
if ( "$count" != "0" ) set tangent_flag = "-DCOST_TANGENT_HAS_NSUB_SLOT"

set adjoint_flag = ""
# count matches in the pipeline; never capture the prototype text itself
# (its trailing brace breaks tcsh word parsing).  Tapenade wraps long
# prototypes, and may break "int *" and "n_sub_settingb" across lines,
# so join the prototype onto one line before testing for the pointer.
set count_ptr = `grep -A12 "void ${head_function}_b *(" $gen_b/tapenade_cost_primal_b.c | sed -n '1,/)/p' | tr '\n' ' ' | grep -c 'int [*] *n_sub_settingb'`
set count_any = `grep -A12 "void ${head_function}_b *(" $gen_b/tapenade_cost_primal_b.c | sed -n '1,/)/p' | grep -c "n_sub_settingb"`
if ( "$count_ptr" != "0" ) then
    set adjoint_flag = "-DCOST_ADJOINT_NSUB_SLOT_BY_POINTER"
else if ( "$count_any" != "0" ) then
    set adjoint_flag = "-DCOST_ADJOINT_NSUB_SLOT_BY_VALUE"
endif
echo "(integer slot flags: tangent '$tangent_flag'  adjoint '$adjoint_flag')"

echo ""
echo "==== 5. compiling"
set objects = ""
foreach gen_c ( $gen_d/*_d.c $gen_b/*_b.c )
    if ( "$gen_c:t" == "tapenade_cost_primal_d.c" ) continue
    if ( "$gen_c:t" == "tapenade_cost_primal_b.c" ) continue
    set obj = $gen_c:r.o
    cc $cflags -I$gen_c:h -I$work_dir -I$kit_dir -c $gen_c -o $obj
    if ( $status != 0 ) then
        echo "ERROR: compiling $gen_c failed"
        exit 1
    endif
    set objects = "$objects $obj"
end

# adStack.c calls clock_gettime(), which strict -std=c11 hides unless a
# POSIX feature level is requested.
cc -std=c11 -O2 -D_POSIX_C_SOURCE=200809L -I$kit_dir -c $kit_dir/adStack.c -o $work_dir/adStack.o
if ( $status != 0 ) then
    echo "ERROR: compiling adStack.c failed"
    exit 1
endif

cc $cflags $tangent_flag -I$gen_d -I$work_dir -I$kit_dir \
    -c tapenade_cost_tangent_adapter.c -o $work_dir/tangent_adapter.o
if ( $status != 0 ) then
    echo "ERROR: compiling the tangent adapter failed"
    exit 1
endif
cc $cflags $adjoint_flag -I$gen_b -I$work_dir -I$kit_dir \
    -c tapenade_cost_adjoint_adapter.c -o $work_dir/adjoint_adapter.o
if ( $status != 0 ) then
    echo "ERROR: compiling the adjoint adapter failed (generated signature differs from expected)"
    exit 1
endif

cc $cflags -Iinclude -c tapenade_cost_primal.c -o $work_dir/tapenade_cost_primal.o
if ( $status != 0 ) exit 1
cc $cflags -Iinclude -c tapenade_adjoint_driver.c -o $work_dir/tapenade_adjoint_driver.o
if ( $status != 0 ) exit 1

cc -o $exe $work_dir/tapenade_adjoint_driver.o $work_dir/tapenade_cost_primal.o \
    $work_dir/tangent_adapter.o $work_dir/adjoint_adapter.o $objects $work_dir/adStack.o \
    build/lib/libsoil.a -lm -Wl,--allow-multiple-definition
if ( $status != 0 ) then
    echo "ERROR: link failed"
    exit 1
endif

echo ""
echo "==== 6. running experiments"
rm -rf $results_dir
mkdir -p $results_dir
./$exe $results_dir forcing/rain_pet_example.csv $experiments | tee $results_dir/results.txt
echo ""
echo "Results: $results_dir/results.txt"
