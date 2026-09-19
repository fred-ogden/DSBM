The NWM Docs (IIRC) say that percolation to gw out of the soil depends on a field capacity threshold.
This code was modified to test and see if that is true, becuase the Noah-MP code does not seem
to apply that threshold.



FC percolation threshold revision

New command-line option:
    --apply-fc-perc-threshold

Semantics:
  Without the option:
    DSBM   - conductivity-based bottom drainage with no theta_fc cutoff.
    NoahMP - native conductivity-based free drainage with no theta_fc cutoff.

  With the option:
    DSBM   - drainage only from water above theta_fc; bottom disc cannot
             be drained below theta_fc.
    NoahMP - same theta_fc threshold and available-water cap applied to
             its bottom free-drainage flux.

This makes it possible to compare:
  no flag : no-FC-threshold versus no-FC-threshold
  flag    : FC-threshold versus FC-threshold

The existing perc_limiter_0_to_1 multiplier remains active in both cases.

Files in this archive are complete revised files, not diffs.

The two pairs of .csh files run dsbm and noah-mp soil moisture routines with/without
this threshold.

Without the fc threshold dsbm and noah-mp soil moisture transport routines in terms
of both states and fluxes agree with each other with NSE and KGE>0.98.
