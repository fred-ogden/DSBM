# DSBM Tapenade soil-parameter and equal-head validation

The two primal wrappers exercise different differentiation paths:

- `tapenade_soil_cost_primal.c` computes a one-year DSBM twin-experiment cost function (normalized squared errors in hourly lateral flow and percolation). Tapenade differentiates it with respect to lateral-flow coefficient `klf_m_per_h`, percolation limiter, saturated conductivity `K_sat_cm_per_h`, and Clapp-Hornberger exponent `b_exp`. The Ksat-to-saturated-capillary-head regression and derived field capacity and wilting point are **inside** the differentiated calculation. Both model and synthetic-observation routines explicitly enable the equal-head limiter.
- `equal_head_ad_primal.c` isolates a downward transfer of 0.015 m from an initially wetter upper disc. The limiter is active, and Tapenade differentiates its Newton solver with respect to `b_exp`. This checks the actual active branch rather than relying on it to occur during a longer simulation.

## One-command validation

From the DSBM repository root:

```tcsh
./tapenade/soil_params/validate_equal_head.csh
```

Requires `cc`, `make`, `python3`, Tapenade 3.16 on `PATH`, and Tapenade's `ADFirstAidKit` at the path configured near the top of the script. The script uses the existing `run_tapenade_soil_params.csh ab` for E10a/E10b, then stages and differentiates the isolated limiter, compiles the tangent and adjoint test drivers, and runs numerical acceptance checks. All generated code is under ignored `tapenade_input/`; outputs are under ignored `output/`.

The script exits **0** only if all checks pass. A nonzero exit indicates a build, generation, execution, or numerical-validation failure. Full-model build logs are in `output/equal_head_validation/full_model_build.log`; E10 results in `output/tapenade_soil_params/results.txt`; active-limiter output in `output/equal_head_validation/{tangent,adjoint}.txt`.

## Interpretation

E10a compares all four cost-function derivatives from tangent and reverse mode at 4 and 12 substeps and tests four directional dot products at each substep count. Acceptance requires relative differences below `1e-9`. E10b compares adjoint gradients with central finite differences across five relative step sizes; at least one step must give relative error below `1e-5` for each parameter. The finite-difference error need not decrease monotonically: large steps have truncation error; small steps can encounter roundoff or branch changes.

The isolated test requires the limited transfer to lie strictly between zero and 0.015 m, an equal-total-head residual below `1e-8` m, a tangent/finite-difference relative error below `1e-7` at one tested step, and tangent/adjoint agreement within `1e-10` relative error.

These tests validate local numerical derivatives at the chosen states and parameters. They do not prove smoothness across limiter activation thresholds, global differentiability, or calibration identifiability from observed streamflow.

## Existing standalone experiments

```tcsh
./tapenade/soil_params/run_tapenade_soil_params.csh ab
./tapenade/soil_params/run_tapenade_soil_params.csh abcd
```

The first runs E10a/E10b derivative checks. The second additionally runs E10c timing/memory and E10d BFGS twin calibration.

FLO 10/26: This validation procedure preserves the verified active-limiter regression point and checks both Tapenade modes against independent finite differences.
