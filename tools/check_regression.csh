#!/bin/tcsh -f
#
# check_regression.csh
#
# Byte-for-byte regression check of the current commit (HEAD) against the
# last accepted baseline, for changes in which output is expected NOT to
# change (renames, refactoring, performance work).  FLO 11/06/2026.
#
# Baseline: 3bcfc9f, infiltration as a flux upper boundary condition (rain
# cascade, control_soil_infiltration_flux_upper_boundary_true_false).  When
# a change that is meant to alter output is accepted, make its commit the
# new baseline_commit below.
#
# Every run script is run in both commits, the same script on both sides:
#   DSBM:    run_dsbm.csh, run_dsbm_analytic.csh, run_dsbm_fc_limit.csh,
#            run_dsbm_infiltration_excess.csh,
#            run_dsbm_analytic_infiltration_excess.csh
#   Noah-MP: run_noahmp.csh, run_noahmp_fc_limit.csh, run_noahmp_fidelity1.csh
# The theta, flux and volume-balance time series must match byte for
# byte; the summary files must match apart from the wall-clock time line.
#
# Both commits are built and run in their own git worktrees under build/
# (ignored by git), so output/ in this checkout is not touched.  Results
# and run logs go to output/regression_check/.
#
# Usage, from the repository root (or from tools/ as ./check_regression.csh):
#     tools/check_regression.csh                (baseline commit 3bcfc9f)
#     tools/check_regression.csh <commit>       (any other baseline that has all the run scripts)
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
set baseline_commit = 3bcfc9f
if ( $#argv > 0 ) set baseline_commit = "$1"
set repo_root    = `pwd`
set baseline_dir = $repo_root/build/regression_baseline
set current_dir  = $repo_root/build/regression_current
set results_dir  = $repo_root/output/regression_check
set runs         = "dsbm dsbm_analytic dsbm_fc_limit dsbm_infiltration_excess dsbm_analytic_infiltration_excess noahmp noahmp_fc_limit noahmp_fidelity1"
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
foreach run ( $runs )
    if ( ! -f $baseline_dir/run_$run.csh ) then
        echo "ERROR: run_$run.csh not in baseline $baseline_commit (predates it?)"
        exit 1
    endif
    if ( ! -f $current_dir/run_$run.csh ) then
        echo "ERROR: run_$run.csh not in HEAD"
        exit 1
    endif
end

echo "==== 2. building both"
foreach worktree_dir ( $baseline_dir $current_dir )
    ( cd $worktree_dir && make -s ) >& $worktree_dir/build.log
    if ( $status != 0 ) then
        echo "ERROR: build failed in $worktree_dir; see $worktree_dir/build.log"
        exit 1
    endif
end

echo "==== 3. running every script in both"
foreach run ( $runs )
    foreach side ( baseline current )
        if ( "$side" == "baseline" ) then
            set worktree_dir = $baseline_dir
        else
            set worktree_dir = $current_dir
        endif
        ( cd $worktree_dir && ./run_$run.csh ) >& $results_dir/${side}_$run.log
        if ( $status != 0 ) then
            echo "ERROR: $side run_$run.csh failed; see $results_dir/${side}_$run.log"
            exit 1
        endif
    end
end

echo "==== 4. comparing"
set n_fail = 0
foreach run ( $runs )
    foreach csv_file ( $csv_files )
        cmp -s $baseline_dir/output/$run/$csv_file $current_dir/output/$run/$csv_file
        if ( $status == 0 ) then
            echo "PASS  $run/$csv_file"
        else
            echo "FAIL  $run/$csv_file"
            @ n_fail = $n_fail + 1
        endif
    end
    grep -v wall-clock $baseline_dir/output/$run/volbal_summary.out >! $results_dir/baseline_${run}_summary.txt
    grep -v wall-clock $current_dir/output/$run/volbal_summary.out >! $results_dir/current_${run}_summary.txt
    cmp -s $results_dir/baseline_${run}_summary.txt $results_dir/current_${run}_summary.txt
    if ( $status == 0 ) then
        echo "PASS  $run/volbal_summary.out (wall-clock line excluded)"
    else
        echo "FAIL  $run/volbal_summary.out (wall-clock line excluded)"
        @ n_fail = $n_fail + 1
    endif
end

if ( $n_fail != 0 ) then
    echo "ERROR: $n_fail comparison(s) failed.  Outputs are in $baseline_dir/output and $current_dir/output"
    exit 1
endif
echo "All comparisons passed.  Outputs are in $baseline_dir/output and $current_dir/output"
echo "Done."
