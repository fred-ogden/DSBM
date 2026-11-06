#!/bin/tcsh -f
#
# check_infiltration_regression.csh
#
# Regression check of the infiltration flux upper boundary condition
# (config keyword control_soil_infiltration_flux_upper_boundary_true_false,
# FLO 11/06/2026) against an earlier commit, where every DSBM run put rain
# into disc 1 only.
#
#   1. FALSE must be bitwise identical to the baseline:
#        baseline run_dsbm.csh          vs current run_dsbm_infiltration_excess.csh
#        baseline run_dsbm_analytic.csh vs current run_dsbm_analytic_infiltration_excess.csh
#      Theta, flux and volume-balance time series byte for byte; summary
#      files apart from the wall-clock line and the new "Infiltration
#      flux upper boundary" line.
#   2. Noah-MP must be bitwise identical (its kernel ignores the keyword):
#        run_noahmp.csh (fidelity 2) and run_noahmp_fidelity1.csh.
#   3. TRUE (current run_dsbm.csh and run_dsbm_analytic.csh) must match
#      the baseline up to the first hour in which the baseline had
#      surface excess (the first hour the disc-1 cap acted); after that
#      the soil moisture differs, as expected.
#   4. A table of totals for every run: rain input, excess, infiltrated
#      rain, AET, lateral flow, percolation, volume balance residual.
#
# Both commits are built and run in their own git worktrees under build/
# (ignored by git), so output/ in this checkout is not touched.  Results
# and run logs go to output/infiltration_regression/.
#
# Usage, from the repository root (or from tools/ as ./check_infiltration_regression.csh):
#     tools/check_infiltration_regression.csh               (baseline 70eb078)
#     tools/check_infiltration_regression.csh <commit>      (any commit before the cascade)
#
# Exit status 0 only if every check passes.  ASCII only.

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
set baseline_commit = 70eb078
if ( $#argv > 0 ) set baseline_commit = "$1"
set repo_root    = `pwd`
set baseline_dir = $repo_root/build/infiltration_baseline
set current_dir  = $repo_root/build/infiltration_current
set results_dir  = $repo_root/output/infiltration_regression
set csv_files    = "theta_timeseries.csv fluxes_timeseries.csv volbal_timeseries.csv"
set baseline_runs = "dsbm dsbm_analytic noahmp noahmp_fidelity1"
set current_runs  = "dsbm_infiltration_excess dsbm_analytic_infiltration_excess dsbm dsbm_analytic noahmp noahmp_fidelity1"
# bitwise pairs "baseline_run current_run"
set bitwise_pairs = ( "dsbm dsbm_infiltration_excess" \
                      "dsbm_analytic dsbm_analytic_infiltration_excess" \
                      "noahmp noahmp" \
                      "noahmp_fidelity1 noahmp_fidelity1" )
set true_runs    = "dsbm dsbm_analytic"
set summary_skip = ( -e wall-clock -e Infiltration )   # "Infiltration flux upper boundary" line (not "Infiltrated")

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
if ( ! -f $current_dir/run_dsbm_infiltration_excess.csh ) then
    echo "ERROR: HEAD has no run_dsbm_infiltration_excess.csh; commit the cascade first"
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

echo "==== 3. running"
foreach run ( $baseline_runs )
    ( cd $baseline_dir && ./run_$run.csh ) >& $results_dir/baseline_$run.log
    if ( $status != 0 ) then
        echo "ERROR: baseline run_$run.csh failed; see $results_dir/baseline_$run.log"
        exit 1
    endif
end
foreach run ( $current_runs )
    ( cd $current_dir && ./run_$run.csh ) >& $results_dir/current_$run.log
    if ( $status != 0 ) then
        echo "ERROR: current run_$run.csh failed; see $results_dir/current_$run.log"
        exit 1
    endif
end

set n_fail = 0

echo ""
echo "==== 4. bitwise checks (FALSE against the baseline; Noah-MP unchanged)"
foreach pair ( $bitwise_pairs:q )
    set names = ( $pair )
    set baseline_out = $baseline_dir/output/$names[1]
    set current_out  = $current_dir/output/$names[2]
    foreach csv_file ( $csv_files )
        cmp -s $baseline_out/$csv_file $current_out/$csv_file
        if ( $status == 0 ) then
            echo "PASS  baseline $names[1] = current $names[2] : $csv_file"
        else
            echo "FAIL  baseline $names[1] = current $names[2] : $csv_file"
            @ n_fail = $n_fail + 1
        endif
    end
    grep -v $summary_skip $baseline_out/volbal_summary.out >! $results_dir/baseline_$names[1]_summary.txt
    grep -v $summary_skip $current_out/volbal_summary.out  >! $results_dir/current_$names[2]_summary.txt
    cmp -s $results_dir/baseline_$names[1]_summary.txt $results_dir/current_$names[2]_summary.txt
    if ( $status == 0 ) then
        echo "PASS  baseline $names[1] = current $names[2] : volbal_summary.out (wall-clock and boundary lines excluded)"
    else
        echo "FAIL  baseline $names[1] = current $names[2] : volbal_summary.out (wall-clock and boundary lines excluded)"
        @ n_fail = $n_fail + 1
    endif
end

echo ""
echo "==== 5. TRUE against the baseline: identical until the disc-1 cap first acted"
foreach run ( $true_runs )
    set baseline_volbal = $baseline_dir/output/$run/volbal_timeseries.csv
    set current_theta   = $current_dir/output/$run/theta_timeseries.csv
    set baseline_theta  = $baseline_dir/output/$run/theta_timeseries.csv
    # line number of the first data line with excess_m > 0 in the baseline
    # (whitespace-delimited; the date and time are fields 1 and 2, excess_m field 4)
    set first_excess_line = `awk 'substr($0, 1, 1) == "#" {next} $4 + 0.0 > 0.0 {print NR; exit}' $baseline_volbal`
    set n_excess_hours    = `awk 'substr($0, 1, 1) == "#" {next} $4 + 0.0 > 0.0 {n++} END {print n + 0}' $baseline_volbal`
    # first line at which the theta time series differ
    set first_diff_line = `cmp $baseline_theta $current_theta |& awk '/differ/ {print $NF; exit}'`
    if ( "$first_diff_line" == "" ) set first_diff_line = none
    echo "$run : baseline hours with excess $n_excess_hours, first at line $first_excess_line; first differing theta line $first_diff_line"
    if ( "$first_diff_line" == "none" ) then
        if ( "$n_excess_hours" != "0" ) then
            echo "FAIL  $run : the baseline had excess, but TRUE output is identical"
            @ n_fail = $n_fail + 1
        else
            echo "PASS  $run : no excess in the baseline, and TRUE output is identical"
        endif
    else if ( "$first_excess_line" == "" ) then
        echo "FAIL  $run : TRUE output differs, but the baseline had no excess"
        @ n_fail = $n_fail + 1
    else if ( $first_diff_line < $first_excess_line ) then
        echo "FAIL  $run : TRUE output differs before the first hour with excess"
        @ n_fail = $n_fail + 1
    else
        echo "PASS  $run : TRUE output identical until the cap first acted"
    endif
end

echo ""
echo "==== 6. totals (m) over the record"
set table = $results_dir/totals.txt
printf "%-44s %10s %10s %10s %10s %10s %10s %12s\n" run rain excess infiltr AET lateral percol residual >! $table
foreach entry ( "baseline dsbm" "baseline dsbm_analytic" "current dsbm_infiltration_excess" "current dsbm_analytic_infiltration_excess" \
                "current dsbm" "current dsbm_analytic" "current noahmp" "current noahmp_fidelity1" )
    set names = ( $entry )
    if ( "$names[1]" == "baseline" ) then
        set summary_file = $baseline_dir/output/$names[2]/volbal_summary.out
    else
        set summary_file = $current_dir/output/$names[2]/volbal_summary.out
    endif
    awk -v label="$names[1] $names[2]" -F= ' \
        /Precipitation input/        {rain = $2 + 0.0} \
        /Surface precipitation excess/ {excess = $2 + 0.0} \
        /Infiltrated precipitation/  {infil = $2 + 0.0} \
        /Evapotranspiration/         {aet = $2 + 0.0} \
        /Lateral flow generated/     {lat = $2 + 0.0} \
        /Percolation to g.w./        {perc = $2 + 0.0} \
        /Volume balance residual/    {res = $2 + 0.0} \
        END {printf "%-44s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %12.4e\n", label, rain, excess, infil, aet, lat, perc, res}' \
        $summary_file >> $table
end
cat $table
echo "(dsbm: lookup table; dsbm_analytic: analytic; _infiltration_excess: FALSE; others TRUE;"
echo " noahmp: fidelity 2, with DSBM lateral flow; noahmp_fidelity1: Noah-MP as written)"

echo ""
if ( $n_fail != 0 ) then
    echo "ERROR: $n_fail check(s) failed.  Outputs are in $baseline_dir/output and $current_dir/output"
    exit 1
endif
echo "All checks passed.  Totals: $table"
echo "Done."
