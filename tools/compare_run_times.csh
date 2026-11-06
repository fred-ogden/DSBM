#!/bin/tcsh -f
#
# compare_run_times.csh
#
# Run-time comparison of the current commit (HEAD) against an earlier
# commit (FLO 11/06/2026).  The default baseline is 0c42e3b, the last
# commit before the equal-head limit on downward transfers was added,
# so the comparison shows the cost of the limiter and of anything since.
#
# Four DSBM configurations are timed, all over the full forcing record
# with configs/soil_params.dat:
#
#   dsbm            run_dsbm.csh: lookup table, full output
#   dsbm_analytic   run_dsbm_analytic.csh: analytic Clapp-Hornberger,
#                   full output
#   kernel_lut      lookup table, --verbosity 0 and only --write-volbal
#   kernel_analytic analytic, --verbosity 0 and only --write-volbal
#
# The full-output runs format several CSV files and one screen line
# per hour, which takes a large share of their time.  The kernel runs
# write only the short volume-balance line per hour (--write-volbal is
# kept because the summary file carrying the wall-clock time is written
# only with it), so they measure mainly the DSBM kernel, the cost that
# matters inside a calibration.
#
# Both commits are built in their own git worktrees under build/
# (ignored by git), so output/ in this checkout is not touched.  Each
# configuration is run n_repeats times per commit, alternating between
# the commits so that machine load affects both alike.  The time
# compared is the driver's "Simulation wall-clock time" from
# volbal_summary.out (the time loop only, not program start-up).  The
# minimum over the repeats is the least noisy measure of cost.
#
# Usage, from the repository root (or from tools/ as ./compare_run_times.csh):
#     tools/compare_run_times.csh                     (baseline 0c42e3b, 5 repeats)
#     tools/compare_run_times.csh <commit>            (any earlier commit)
#     tools/compare_run_times.csh <commit> <repeats>
#
# Output: output/run_times/run_times.txt.  ASCII only.

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
set baseline_commit = 0c42e3b
set n_repeats = 5
if ( $#argv > 0 ) set baseline_commit = "$1"
if ( $#argv > 1 ) set n_repeats = "$2"
set repo_root    = `pwd`
set baseline_dir = $repo_root/build/timing_baseline
set current_dir  = $repo_root/build/timing_current
set results_dir  = $repo_root/output/run_times
set runs         = "dsbm dsbm_analytic kernel_lut kernel_analytic"
# options shared by the kernel runs (paths relative to each worktree)
set kernel_options = "--config configs/soil_params.dat --solver dsbm --forcing forcing/rain_pet_example.csv --verbosity 0 --write-volbal --timestamp datetime --dz 0.1,0.3,0.6,1.0 --zc 0.05,0.25,0.7,1.5"
# lookup-table options: commits before FLO 11/06/2026 took a table size and
# minimum (400 points from 1e-6 were used); later commits have a fixed
# 2-point table and reject those options.  Chosen per worktree below.
set lut_options_old = "--use-lut --lut-n 400 --lut-Theta-min 1e-6"
set lut_options_new = "--use-lut"

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

echo "==== 3. timing, $n_repeats repeats per configuration per commit"
set report = $results_dir/run_times.txt
echo "Run times (s), driver wall-clock of the time loop; baseline $baseline_commit, current HEAD" >! $report
foreach run ( $runs )
    echo "---- $run"
    set times_baseline = ( )
    set times_current = ( )
    set repeat = 0
    while ( $repeat < $n_repeats )
        @ repeat = $repeat + 1
        foreach worktree_dir ( $baseline_dir $current_dir )
            if ( "$run" == "kernel_lut" ) then
                set outdir = output/timing_$run
                mkdir -p $worktree_dir/$outdir
                set lut_options = "$lut_options_new"
                grep -q "&o->lut_n" $worktree_dir/src/bmi_soil_driver.c
                if ( $status == 0 ) set lut_options = "$lut_options_old"
                ( cd $worktree_dir && bin/soil_driver $kernel_options $lut_options --outdir $outdir ) >& /dev/null
                set run_status = $status
            else if ( "$run" == "kernel_analytic" ) then
                set outdir = output/timing_$run
                mkdir -p $worktree_dir/$outdir
                ( cd $worktree_dir && bin/soil_driver $kernel_options --outdir $outdir ) >& /dev/null
                set run_status = $status
            else
                set outdir = output/$run
                ( cd $worktree_dir && ./run_$run.csh ) >& /dev/null
                set run_status = $status
            endif
            if ( $run_status != 0 ) then
                echo "ERROR: $run failed in $worktree_dir"
                exit 1
            endif
            set seconds = `grep "wall-clock" $worktree_dir/$outdir/volbal_summary.out | awk '{print $(NF-1)}'`
            if ( "$worktree_dir" == "$baseline_dir" ) then
                set times_baseline = ( $times_baseline $seconds )
            else
                set times_current = ( $times_current $seconds )
            endif
        end
    end
    set min_baseline = `echo $times_baseline | awk '{m=$1; for(i=2;i<=NF;i++) if($i<m) m=$i; print m}'`
    set min_current  = `echo $times_current  | awk '{m=$1; for(i=2;i<=NF;i++) if($i<m) m=$i; print m}'`
    set ratio = `echo $min_baseline $min_current | awk '{if ($1 > 0) printf "%.3f", $2/$1; else print "n/a"}'`
    echo "" >> $report
    echo "$run" >> $report
    echo "   baseline : $times_baseline   minimum $min_baseline" >> $report
    echo "   current  : $times_current   minimum $min_current" >> $report
    echo "   current/baseline (minima) = $ratio" >> $report
end

cat $report
echo "Report written to $report"
echo "Done."
