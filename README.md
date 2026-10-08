# DSBM: Discretized Soil-Moisture Balance Model

DSBM is a compact, stateless soil-moisture model that simulates vertical
soil-moisture dynamics in a **vertically discretized soil column having
uniform soil hydraulic characteristics**.

The soil column is divided numerically into a small number of
finite-thickness **discs**. These discs are computational control
volumes, not soil layers representing different soil types or horizons.
The same soil hydraulic parameters apply throughout the column.

This repository also contains an independent C implementation of the
Noah-MP soil-moisture solver for direct numerical comparison with DSBM.
The comparison implementation follows the Noah-MP
`SOILWATER -> SRT -> SSTEP -> ROSR12` solution path and provides a
reference against which the simpler DSBM formulation can be evaluated.

The immediate purposes of this repository are to:

1.  preserve and develop the DSBM formulation;
2.  provide a reproducible comparison of DSBM with the Noah-MP
    soil-moisture solution using common soil hydraulic properties,
    soil-column geometry, initial soil-moisture storage, and water
    forcing; and
3.  stage the differentiation of DSBM with the Tapenade automatic
    differentiation tool, as a step toward gradient-based calibration of
    CFE3.1 (see [Differentiability](#differentiability-tapenade-experiments)).

## DSBM formulation

DSBM represents a soil column of uniform hydraulic characteristics by
dividing it vertically into discrete control volumes, referred to here
as **discs**.

The standard four-disc configuration used in the supplied comparison is
consistent with the discretization used in Noah-MP as applied in the
National Water Model version 3.1 and earlier:

| Disc | Thickness (m) | Center depth (m) |
|-----:|--------------:|-----------------:|
|    1 |          0.10 |             0.05 |
|    2 |          0.30 |             0.25 |
|    3 |          0.60 |             0.70 |
|    4 |          1.00 |             1.50 |

The total soil depth is 2.0 m.

The different disc thicknesses provide vertical numerical
discretization. They do **not** imply changes in soil texture, soil
type, or hydraulic properties with depth. A single set of soil hydraulic
parameters describes the entire soil column.

Soil moisture is represented by volumetric soil moisture content,
theta, in each disc. Soil moisture is redistributed between adjacent
discs according to Darcy-Buckingham fluxes. Positive vertical flux is
downward.

## Soil hydraulic relationships

DSBM uses Clapp-Hornberger relationships to calculate matric potential
and unsaturated hydraulic conductivity from volumetric soil moisture.

For a uniform soil,

```math
\psi(\theta) = \psi_{sat} \left( \frac{\theta}{\theta_{sat}} \right)^{-b}
```

and

```math
K(\theta) = K_{sat} \left( \frac{\theta}{\theta_{sat}} \right)^{2b+3} .
```

The parameters theta_sat, K_sat, psi_sat and b characterize the soil
and are common to all discs.

For an interface between adjacent discs, DSBM evaluates the
Darcy-Buckingham flux using the hydraulic states of the two adjoining
control volumes and their separation distance. Thus, spatial variability
in theta, psi, and K develops dynamically even though the underlying
soil hydraulic parameterization is uniform with depth.

Each hourly timestep is divided into substeps (4 by default; see
[Substeps per timestep](#substeps-per-timestep)). Within a substep the
transfers between discs are capped by the soil moisture available in the
donor disc and the storage available in the receiving discs, which keeps
the solution stable and mass-conserving.

An optional lookup table can be used for the Clapp-Hornberger hydraulic
functions.

## Noah-MP comparison solver

The repository includes `src/noahmp_soilmoisture_stateless.c`, an
independent C implementation of the Noah-MP vertical soil-moisture
solution used for comparison with DSBM.

For the comparison experiments in this repository, Noah-MP is configured
with the same uniform soil hydraulic characteristics throughout the soil
column and the same vertical discretization used by DSBM. The objective
is to compare the numerical soil-moisture formulations rather than to
compare different soil profiles.

The Noah-MP implementation retains the essential numerical structure of
the original solver:

-   calculation of soil hydraulic diffusivity and conductivity;
-   assembly of the implicit soil-moisture equations;
-   tridiagonal matrix solution;
-   update of disc soil moisture; and
-   conductivity-based lower-boundary drainage.

The two formulations can therefore be run with the same forcing, soil
hydraulic parameters, discretization, and initial soil-moisture storage.

Select the formulation with:

``` text
--solver dsbm
```

or

``` text
--solver noahmp
```

## Lower-boundary percolation and field capacity

A central purpose of the comparison is to distinguish differences caused
by the numerical treatment of vertical soil-moisture movement from
differences caused by the lower-boundary condition.

By default, both comparison implementations use conductivity-based
bottom drainage without imposing field capacity as a percolation
threshold.

The optional switch

``` text
--apply-fc-perc-threshold
```

changes the lower-boundary treatment so that percolation occurs only
from soil moisture above field capacity and removal of percolating water
does not reduce the bottom-disc moisture below `theta_fc`.

For Noah-MP, this switch is an experimental comparison option. It is
**not** the native free-drainage behavior represented by the Noah-MP
`OPT_RUN = 3` and `OPT_RUN = 7` soil-moisture paths associated with the
Schaake and Xinanjiang surface-water partitioning configurations used by
the National Water Model.

The option is retained because it permits controlled comparisons of DSBM
and Noah-MP under either lower-boundary assumption. See also
`README_FC_THRESHOLD.txt`.

## Initial DSBM and Noah-MP comparison

For the current test case, when the field-capacity percolation threshold
is **not** applied, DSBM and the Noah-MP soil-moisture solution agree
very closely. NSE and KGE values for the simulated soil-moisture states
and water fluxes are generally greater than 0.98.

This result is notable because the two formulations solve the vertical
soil-moisture problem using substantially different numerical methods:

-   DSBM explicitly transfers soil moisture between finite soil discs
    using Darcy-Buckingham fluxes and sub-timesteps.
-   Noah-MP solves its discretized soil-moisture equations using an
    implicit tridiagonal matrix solution.

Both calculations in this experiment represent the same vertically
uniform soil hydraulic properties and use comparable boundary
conditions.

The reported NSE and KGE values describe the supplied comparison
experiment. They should not be interpreted as validation across
arbitrary soil properties, forcings, initial conditions, climates, or
numerical configurations.

The comparison figures and summaries in `results/` were produced before
the current defaults (fixed 4 DSBM substeps, exact exponential lateral
removal, computed wilting point) were adopted. They can be reproduced exactly; see
[Reproducing the original DSBM output](#reproducing-the-original-dsbm-output).

## Repository layout

``` text
.
|-- Makefile
|-- README.md, README.1st, README_FC_THRESHOLD.txt, LICENSE
|-- configs/
|   `-- soil_params.dat                  soil parameters (key=value)
|-- forcing/
|   `-- rain_pet_example.csv             hourly rainfall and PET, 78168 hours
|-- include/
|   |-- dsbm_census_hooks.h              branch-census hooks (empty in normal builds)
|   |-- dsbm_soilmoisture_stateless.h
|   |-- noahmp_soilmoisture_stateless.h
|   |-- soil_cli.h
|   |-- soil_config.h
|   |-- soil_data_types.h
|   |-- soil_helpers.h
|   `-- util.h
|-- src/
|   |-- bmi_soil_driver.c                command-line driver
|   |-- dsbm_census.c                    branch-census counters (census build only)
|   |-- dsbm_soilmoisture_stateless.c    DSBM kernel
|   |-- noahmp_soilmoisture_stateless.c  Noah-MP comparison kernel
|   |-- soil_helpers.c                   Clapp-Hornberger, fluxes, lateral removal
|   `-- util.c
|-- results/                             DSBM vs Noah-MP comparison figures
|-- run_dsbm.csh                         DSBM, lookup table
|-- run_dsbm_analytic.csh                DSBM, analytic Clapp-Hornberger
|-- run_dsbm_fc_limit.csh                DSBM, field-capacity percolation threshold
|-- run_noahmp.csh                       Noah-MP
|-- run_noahmp_fc_limit.csh              Noah-MP, field-capacity percolation threshold
|-- tools/
|   |-- run_branch_census.csh            branch-activity census
|   `-- run_nsub_timing.csh              substep count versus cost and results
`-- tapenade/                            automatic differentiation experiments
    |-- initial/                         first one-timestep experiment
    |-- klf/                             E1-E7: lateral rate constant k_lf
    |-- cost/                            E8-E9: cost-function gradient, adjoint
    `-- soil_params/                     E10: four parameters, phi_sat from Ksat
```

Generated executables, build products, model output (`output/`), and
Tapenade working files (`tapenade_input/`) are intentionally excluded
from version control.

## Building

Build the model with:

``` bash
make
```

The executable is written to:

``` text
bin/soil_driver
```

The build also creates the static soil library:

``` text
build/lib/libsoil.a
```

For a debug build:

``` bash
make debug
```

To remove generated build products:

``` bash
make veryclean
```

The number of soil discs and minimum soil moisture can be overridden at
build time if required:

``` bash
make NDISC=4 THETA_MIN=1.0e-03
```

Changing `NDISC` changes the numerical discretization of the soil
column; it does not define different soil hydraulic properties for the
individual discs.

A separate census build (see [Branch-activity census](#branch-activity-census))
is made with:

``` bash
make census
```

## Running the comparison

Five convenience scripts are supplied:

``` text
run_dsbm.csh
run_dsbm_analytic.csh
run_noahmp.csh
run_dsbm_fc_limit.csh
run_noahmp_fc_limit.csh
```

The first group runs the conductivity-based lower-boundary comparison
(`run_dsbm.csh` uses the Clapp-Hornberger lookup table,
`run_dsbm_analytic.csh` the analytic functions):

``` bash
./run_dsbm.csh
./run_dsbm_analytic.csh
./run_noahmp.csh
```

The second pair applies the optional field-capacity percolation
threshold to both formulations:

``` bash
./run_dsbm_fc_limit.csh
./run_noahmp_fc_limit.csh
```

The scripts use the same uniform soil hydraulic properties, vertical
discretization, forcing, and initial total soil-moisture storage so that
state and flux time series can be compared directly. Output goes to
`output/<script name without run_>/`.

The DSBM scripts use the current DSBM defaults. To reproduce the output
from before those defaults changed, see the next section.

## Reproducing the original DSBM output

Three defaults have changed since commit `1a8d1e2`: the DSBM substep
count (fixed at 4, previously adaptive), the DSBM lateral removal (exact
exponential, previously forward Euler), and the wilting point (computed
at 15 atmospheres, previously fixed at 0.10, for both DSBM and Noah-MP).
The original behavior is still available.  Add this line to the config
file:

``` text
soil_wilting_point_m3_per_m3=0.10[V V-1]
```

and, for DSBM runs, these options:

``` text
--nsub 0 --lateral-forward-euler
```

The output then matches the code at commit `1a8d1e2` byte for byte:
theta, flux, and volume-balance time series, and the summary file apart
from its wall-clock time line.  The screen output also matches apart
from one added first line reporting the wilting point.  This has been
checked for DSBM with the lookup table and with the analytic functions,
each with and without `--apply-fc-perc-threshold`, and for Noah-MP with
and without it.

## Driver examples

A DSBM run can be made with:

``` bash
bin/soil_driver \
    --config configs/soil_params.dat \
    --solver dsbm \
    --forcing forcing/rain_pet_example.csv \
    --outdir output/dsbm_example \
    --write-theta \
    --write-fluxes \
    --write-volbal \
    --dz 0.1,0.3,0.6,1.0 \
    --zc 0.05,0.25,0.7,1.5
```

For the corresponding Noah-MP calculation, use `--solver noahmp`.

DSBM additionally supports a Clapp-Hornberger lookup table:

``` text
--use-lut --lut-n 400 --lut-Theta-min 1e-6
```

The initial state is a hydrostatic profile. By default its total storage
comes from `state_soil_reservoir_init_storage_m` in the config file; it
can instead be set with `--init-target-storage <m>` or from a water
table depth with `--init-wt-depth <m>`.

`bin/soil_driver` with an unknown option (for example `--help`) prints
the full option list.

## Forcing format

The driver reads an hourly comma-separated file with a required header
row:

``` text
datetime,rainfall_mm_per_h,potential_et_mm_per_h
2012/10/01 00:00:00,0.000000,0.000000
```

Lines beginning with `#` are ignored.

Rainfall and potential evapotranspiration are supplied in mm per hour.

## Output

The driver can write:

-   volumetric soil-moisture time series for each disc
    (`theta_timeseries.csv`);
-   water-flux time series (`fluxes_timeseries.csv`);
-   timestep volume-balance diagnostics (`volbal_timeseries.csv`); and
-   a volume-balance summary (`volbal_summary.out`), which also reports
    the total number of substeps and the wall-clock time.

Enable these with:

``` text
--write-theta
--write-fluxes
--write-volbal
```

The output time coordinate can be an integer timestep index, Julian
Date, or calendar time (`--timestamp timestep|juliandate|datetime`).

## Soil hydraulic parameters

The driver reads a simple parameter file of key=value pairs, with units
in brackets. `configs/soil_params.dat`:

``` text
soil_depth_m=2.0[m]
soil_Clapp_Hornberger_exponent_b=4.05[]
soil_sat_hydraulic_conductivity_cm_per_h=1.2168[cm h-1]
soil_saturated_capillary_head_calc_from_ksat=TRUE // Calculates Hc(cm) = 10.415*Ksat(cm/h)^(-0.3266)
soil_effective_porosity=0.439[V V-1]
soil_field_capacity_Pcap_over_Patm_0_1=0.333[P P-1]
soil_reservoir_rate_const_to_subsurface_lateral_flow=0.000005[m h-1]
soil_to_gw_percolation_rate_limiter_0_to_1=0.81[]
state_soil_reservoir_init_storage_m=0.7843498367003[m]
```

Some quantities are derived rather than read:

-   With `soil_saturated_capillary_head_calc_from_ksat=TRUE`, the
    saturated capillary head is computed from Ksat with the CFE3.1
    regression `phi_sat (cm) = 10.415 * Ksat (cm/h)^(-0.3266)`, so the
    two inversely related soil properties cannot be set independently.
    Otherwise it is read from `soil_sat_capillary_head_cm`.
-   Field capacity `theta_fc` is the Clapp-Hornberger soil moisture at
    the field-capacity pressure head
    (`soil_field_capacity_Pcap_over_Patm_0_1` times atmospheric pressure
    head), and `theta_aet_eq_pet` (the soil moisture at and above which
    AET equals the assigned PET) is set equal to it.
-   The wilting point `theta_wp` is the Clapp-Hornberger soil moisture
    at a capillary pressure of 15 atmospheres, as in CFE3.1:
    `theta_wp = theta_sat * (15 * psi_atm / phi_sat)^(-1/b)`, which is
    0.0712 for `configs/soil_params.dat`.  The driver prints it at the
    start of a run.  Two optional keys change it:
    `soil_wilting_point_Pcap_over_Patm=<ratio>[P P-1]` changes the
    pressure ratio, and `soil_wilting_point_m3_per_m3=<value>[V V-1]`
    sets theta_wp directly, overriding the calculation.
-   The residual soil moisture `theta_r` is fixed at 0 in the driver.

These are properties or parameters of the **uniform soil column**. The
current DSBM formulation does not assign an independent set of hydraulic
properties to each disc.

## Lateral subsurface flow

Each disc above field capacity loses soil moisture to the lateral
(Nash cascade) reservoir as a linear reservoir:

``` text
d(theta_i)/dt = -klf * (theta_i - theta_fc) / (dz_i * (theta_sat - theta_fc))
```

with no lateral flow from a disc at or below theta_fc.

By default this equation is integrated exactly over each substep:

``` text
theta_i - theta_fc  <-  (theta_i - theta_fc) * exp(-klf * dt_sub / (dz_i * (theta_sat - theta_fc)))
```

Lateral removal by itself can never take a disc below theta_fc, so the
exact solution needs no storage cap, and its derivative with respect to
klf is never artificially zero, which matters for gradient-based
calibration.  Other processes (Darcy-Buckingham drainage, root water
uptake, bare-soil evaporation from disc 1, and percolation from disc 4)
can and do take discs below theta_fc; lateral flow from such a disc is
then zero.  The volume removed by the lateral step equals the change in
storage it causes.

The legacy forward-Euler step, capped at the soil moisture available
above theta_fc, can be selected with:

``` text
--lateral-forward-euler
```

## Substeps per timestep

By default DSBM uses a fixed 4 substeps per hourly timestep.  A fixed
count keeps outputs, and their derivatives with respect to calibration
parameters, free of the jumps that occur when a parameter change moves a
timestep across one of the adaptive scheme's severity thresholds.  The
count can be changed with:

``` text
--nsub N      fixed at N substeps per timestep (N > 0)
--nsub 0      original adaptive choice (1 to 12 substeps per timestep,
              from the Darcy-Buckingham fluxes and rain rate)
```

With `--nsub 0`, `--n-sub-minimum N` sets a floor on the adaptive count,
and `--substep-lateral-severity` adds a lateral-removal criterion to the
adaptive choice.  The total number of substeps and the wall-clock time
are reported in `volbal_summary.out`.

`tools/run_nsub_timing.csh` runs the full record (configs/soil_params.dat,
78168 hours) for a list of settings.  Totals over the record:

| `--nsub`   | substeps per hour | lateral (m) | percolation (m) | AET (m) |
|------------|------------------:|------------:|----------------:|--------:|
| 0 adaptive |              1.18 |      0.5372 |          6.3367 |  4.0149 |
| 1          |              1    |      0.5340 |          6.2145 |  3.9968 |
| 2          |              2    |      0.5373 |          6.3184 |  4.0149 |
| 4 default  |              4    |      0.5384 |          6.3316 |  4.0149 |
| 12         |             12    |      0.5389 |          6.3341 |  4.0149 |
| 24         |             24    |      0.5391 |          6.3343 |  4.0149 |

With 4 substeps, percolation is within 0.05% and lateral flow within
0.15% of the 24-substep values.  A single substep per hour is not
recommended: on this record it gives about 2% less percolation and
0.5% less AET than the finer settings.  Run time grows much more slowly
than the substep count (roughly 0.2 s for 1 substep per hour, 0.3 s for
4, and 0.8 s for 24, for the whole record on one workstation; timings
vary from run to run).

## Branch-activity census

`make census` builds `bin/soil_driver_census`, the same driver and
kernel with counting hooks compiled in.  It reports how often each
threshold, cap, and branch in the DSBM kernel is active over a run:
Darcy-Buckingham flux direction and transfer caps, percolation cap, ET
stress regimes, lateral activity, theta_fc and theta_wp crossings, and
roundoff clamps.  `tools/run_branch_census.csh` runs it for a list of
lateral rate constants:

``` bash
make census
tools/run_branch_census.csh 1e-3 1e-2
```

In the normal build the hooks expand to nothing, so the production
kernel is unchanged.

## Differentiability (Tapenade experiments)

DSBM is being made differentiable with the Tapenade source-to-source
automatic differentiation tool, as a staging step for CFE3.1.  The goal
is the gradient of a calibration cost function with respect to a vector
of soil parameters.  The kernel source has been arranged so that
Tapenade can differentiate it in tangent (forward) and reverse (adjoint)
modes: all locals declared at the top of each function and no ternary
operators.  These changes do not alter the arithmetic.

Each experiment directory under `tapenade/` contains a primal wrapper,
adapters for the generated code, a driver, and a run script.  A run
script can be started from the repository root or from its own
directory; it rebuilds the library, runs Tapenade, compiles, and writes
results to `output/tapenade_<name>/`.  Tapenade itself is required only
for these experiments, not to build or run the model.

| Directory               | Experiments | Content |
|-------------------------|-------------|---------|
| `tapenade/initial/`     | E0          | first one-timestep test: percolation with respect to Ksat |
| `tapenade/klf/`         | E1-E7       | tangent of lateral flow with respect to k_lf; substep and lateral-scheme effects |
| `tapenade/cost/`        | E8-E9       | gradient of a cost function (sum of 1 - NSE for lateral flow and percolation) with respect to k_lf and perc_limiter; tangent and adjoint; BFGS calibration |
| `tapenade/soil_params/` | E10         | four parameters (k_lf, perc_limiter, Ksat, b), with phi_sat computed from Ksat and theta_fc from phi_sat and b inside the differentiated code |

For example:

``` bash
tapenade/soil_params/run_tapenade_soil_params.csh
```

The calibration experiments are twin experiments: the "observations"
are DSBM's own outputs at known parameter values.  They test the
derivative machinery (tangent and adjoint agree to about 1e-14, and BFGS
recovers the true parameters), not whether the parameters can be
identified from observed discharge, which requires CFE3.1 and routing.

## Numerical comparison objective

The DSBM and Noah-MP experiment is designed primarily to answer:

> How closely can a simple discrete soil-moisture balance formulation
> reproduce the vertical soil-moisture states and water fluxes produced
> by the Noah-MP implicit soil-moisture solver when both represent the
> same uniform soil?

This requires keeping the physical problem as nearly identical as
possible between the two formulations:

-   identical uniform soil hydraulic parameters;
-   identical soil-column depth;
-   identical vertical discretization;
-   identical initial soil-moisture storage and moisture profile;
-   identical water forcing; and
-   comparable upper and lower boundary conditions.

Differences in the resulting state and flux time series can then be
attributed primarily to the numerical formulations and to any explicitly
selected process differences.

The optional field-capacity percolation threshold provides one
controlled experiment for separating a lower-boundary assumption from
the behavior of the vertical soil-moisture solution itself.

## Volume conservation

Volume conservation is a primary diagnostic.

Each timestep accounts for rainfall entering the soil, surface excess,
evapotranspiration, lateral flow, bottom percolation, and change in
total soil-moisture storage.  The residual over the supplied 78168-hour
record is of order 1e-13 m.

Comparison experiments should examine the volume-balance residual in
addition to statistical agreement between DSBM and Noah-MP state and
flux time series.

## Status

DSBM is research and development software.

The repository currently contains:

-   the stateless DSBM soil-moisture formulation;
-   an independent C implementation of the Noah-MP soil-moisture
    solution for numerical comparison;
-   a common command-line driver;
-   common forcing and configuration capability;
-   volume-balance diagnostics;
-   reproducible DSBM and Noah-MP comparison scripts;
-   a branch-activity census build; and
-   Tapenade tangent and adjoint experiments on the DSBM kernel.

The initial comparison demonstrates close agreement for the supplied
experiment. Broader evaluation should include different uniform soil
hydraulic properties, initial moisture states, precipitation regimes,
evapotranspiration demand, and numerical timestep conditions.

## Author

Fred L. Ogden\
NOAA/National Weather Service\
Office of Water Prediction

DSBM originated during development of the CFE soil-moisture formulation
and subsequent investigation of a computationally simple, discretized
representation of vertical soil-moisture dynamics.

## License

Source files originating from CFE/DSBM development contain Apache
License 2.0 notices.

Before assigning a repository-wide license, the provenance and
applicable licensing terms of the Noah-MP comparison implementation
should be documented explicitly.
