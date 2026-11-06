#!/bin/tcsh -f
#
# check_regression.csh
#
# Byte-for-byte regression check of the current commit against an
# earlier commit, for the changes in which output is expected NOT to
# change (FLO 11/06/2026):
#
#   DSBM: run_dsbm.csh, run_dsbm_analytic.csh, run_dsbm_fc_limit.csh.
#     Removing the unlimited downward-transfer option leaves the default
#     (equal-head limit) unchanged.
#
#   Noah-MP: the old run_noahmp.csh and run_noahmp_fc_limit.csh, which
#     always applied the forward-Euler lateral sink, against the current
#     driver at --noahmp-fidelity 2 --lateral-forward-euler (same
#     options otherwise).
#
# The theta, flux and volume-balance time series must match byte for
# byte.  The summary files must match apart from the wall-clock time
# line and the Noah-MP fidelity line (new).
#
# Both commits are built and run in their own git worktrees under
# build/ (ignored by git), so output/ in this checkout is not touched.
# Results and run logs go to output/regression_check/.
#
# Usage, from the repository root (or from tools/ as ./check_regression.csh):
#     tools/check_regression.csh                (baseline commit da787c6)
#     tools/check_regression.csh <commit>       (any earlier commit)
#
# FLO 11/06/2026: since the rain cascade (infiltration as a flux upper
# boundary condition) the current run_dsbm*.csh output differs from any
# earlier commit, so this check applies only to commits before it.  For
# the cascade, use tools/check_infiltration_regression.csh.
#
# Exit status 0 only if every comparison passes.  ASCII only.

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
set baseline_commit = da787c6
if ( $#argv > 0 ) set baseline_commit = "$1"
set repo_root    = `pwd`
set baseline_dir = $repo_root/build/regression_baseline
set current_dir  = $repo_root/build/regression_current
set results_dir  = $repo_root/output/regression_check
set dsbm_runs    = "dsbm dsbm_analytic dsbm_fc_limit"
set csv_files    = "theta_timeseries.csv fluxes_timeseries.csv volbal_timeseries.csv"

git rev-parse --verify --quiet "$baseline_commit^{commit}" >& /dev/null
if ( $status != 0 ) then
    echo "ERROR: baseline commit $baseline_commit not found"
    exit 1
endif

rm -rf $results_dir
mkdir -p $results_dir

echo "==== 1. worktrees: baseline $baseline_commit and current HEAD"
foreach worktree_dir ( $baseline_dir $current_dir )
    if ( -d $worktree_dir ) git worktree remove --force $worktree_dir
end
git worktree prune
git worktree add --detach $baseline_dir $baseline_commit >& $results_dir/worktree.log
if ( $status != 0 ) then
    echo "ERROR: could not create baseline worktree; see $results_dir/worktree.log"
    exit 1
endif
git worktree add --detach $current_dir HEAD >>& $results_dir/worktree.log
if ( $status != 0 ) then
    echo "ERROR: could not create current worktree; see $results_dir/worktree.log"
    exit 1
endif

echo "==== 2. building both"
foreach worktree_dir ( $baseline_dir $current_dir )
    ( cd $worktree_dir && make -s ) >& $worktree_dir/build.log
    if ( $status != 0 ) then
        echo "ERROR: build failed in $worktree_dir; see $worktree_dir/build.log"
        exit 1
    endif
end

echo "==== 3. running the baseline scripts"
foreach run ( $dsbm_runs noahmp noahmp_fc_limit )
    ( cd $baseline_dir && ./run_$run.csh ) >& $results_dir/baseline_$run.log
    if ( $status != 0 ) then
        echo "ERROR: baseline run_$run.csh failed; see $results_dir/baseline_$run.log"
        exit 1
    endif
end

echo "==== 4. running the current DSBM scripts and Noah-MP fidelity 2 forward Euler"
foreach run ( $dsbm_runs )
    ( cd $current_dir && ./run_$run.csh ) >& $results_dir/current_$run.log
    if ( $status != 0 ) then
        echo "ERROR: current run_$run.csh failed; see $results_dir/current_$run.log"
        exit 1
    endif
end
# Same options as the old run_noahmp.csh and run_noahmp_fc_limit.csh,
# plus --noahmp-fidelity 2 --lateral-forward-euler.
foreach run ( noahmp noahmp_fc_limit )
    set fc_option = ""
    if ( "$run" == "noahmp_fc_limit" ) set fc_option = "--apply-fc-perc-threshold"
    mkdir -p $current_dir/output/$run
    ( cd $current_dir && bin/soil_driver --config configs/soil_params.dat --solver noahmp --noahmp-fidelity 2 --lateral-forward-euler $fc_option --forcing forcing/rain_pet_example.csv --outdir output/$run --verbosity 1 --write-theta --write-fluxes --write-volbal --timestamp datetime --dz 0.1,0.3,0.6,1.0 --zc 0.05,0.25,0.7,1.5 ) >& $results_dir/current_$run.log
    if ( $status != 0 ) then
        echo "ERROR: current Noah-MP $run run failed; see $results_dir/current_$run.log"
        exit 1
    endif
end

echo "==== 5. comparing"
set n_fail = 0
foreach run ( $dsbm_runs noahmp noahmp_fc_limit )
    foreach csv_file ( $csv_files )
        cmp -s $baseline_dir/output/$run/$csv_file $current_dir/output/$run/$csv_file
        if ( $status == 0 ) then
            echo "PASS  $run/$csv_file"
        else
            echo "FAIL  $run/$csv_file"
            @ n_fail = $n_fail + 1
        endif
    end
    grep -v -e "wall-clock" -e "Noah-MP fidelity" $baseline_dir/output/$run/volbal_summary.out >! $results_dir/baseline_${run}_summary.txt
    grep -v -e "wall-clock" -e "Noah-MP fidelity" $current_dir/output/$run/volbal_summary.out >! $results_dir/current_${run}_summary.txt
    cmp -s $results_dir/baseline_${run}_summary.txt $results_dir/current_${run}_summary.txt
    if ( $status == 0 ) then
        echo "PASS  $run/volbal_summary.out (wall-clock and fidelity lines excluded)"
    else
        echo "FAIL  $run/volbal_summary.out (wall-clock and fidelity lines excluded)"
        @ n_fail = $n_fail + 1
    endif
end

if ( $n_fail != 0 ) then
    echo "ERROR: $n_fail comparison(s) failed.  Outputs are in $baseline_dir/output and $current_dir/output"
    exit 1
endif
echo "All comparisons passed.  Outputs are in $baseline_dir/output and $current_dir/output"
echo "Done."
