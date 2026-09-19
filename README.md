# DSBM --- Discretized Soil-Moisture Balance Model

DSBM is a compact, stateless soil-water model that simulates vertical
soil-moisture dynamics in a **vertically discretized soil column having
uniform soil hydraulic characteristics**.

The soil column is divided numerically into a small number of
finite-thickness **discs**. These discs are computational control
volumes, not soil layers representing different soil types or horizons.
The same soil hydraulic parameters apply throughout the column.

This repository also contains an independent C implementation of the
Noah-MP soil-water solver for direct numerical comparison with DSBM. The
comparison implementation follows the Noah-MP
`SOILWATER -> SRT -> SSTEP -> ROSR12` solution path and provides a
reference against which the simpler DSBM formulation can be evaluated.

The immediate purposes of this repository are to:

1.  preserve and develop the DSBM formulation; and
2.  provide a reproducible comparison of DSBM with the Noah-MP
    soil-water solution using common soil hydraulic properties,
    soil-column geometry, initial soil-water storage, and water forcing.

## DSBM formulation

DSBM represents a soil column of uniform hydraulic characteristics by
dividing it vertically into discrete control volumes, referred to here
as **discs**.

For the standard four-disc configuration used in the supplied
comparison is consistent with those used in Noah-MP as applied in National Water Model version 3.1 and earlier:

  Disc     Thickness (m)   Center depth (m)
  ------ --------------- ------------------
  1                 0.10               0.05
  2                 0.30               0.25
  3                 0.60               0.70
  4                 1.00               1.50

The total soil depth is 2.0 m.

The different disc thicknesses provide vertical numerical
discretization. They do **not** imply changes in soil texture, soil
type, or hydraulic properties with depth. A single set of soil hydraulic
parameters describes the entire soil column.

Soil moisture is represented by volumetric water content,
(`\theta`{=tex}), in each disc. Water is redistributed between adjacent
discs according to Darcy-Buckingham fluxes. Positive vertical flux is
downward.

## Soil hydraulic relationships

DSBM uses Clapp-Hornberger relationships to calculate matric potential
and unsaturated hydraulic conductivity from volumetric soil moisture.

For a uniform soil,

\[ `\psi`{=tex}(`\theta`{=tex}) = `\psi`{=tex}\_{sat}
`\left`{=tex}(`\frac{\theta}{\theta_{sat}}`{=tex}`\right`{=tex})\^{-b}
\]

and

\[ K(`\theta`{=tex}) = K\_{sat}
`\left`{=tex}(`\frac{\theta}{\theta_{sat}}`{=tex}`\right`{=tex})\^{2b+3}.
\]

The parameters (`\theta`{=tex}*{sat}), (K*{sat}), (`\psi`{=tex}\_{sat}),
and (b) characterize the soil and are common to all discs.

For an interface between adjacent discs, DSBM evaluates the
Darcy-Buckingham flux using the hydraulic states of the two adjoining
control volumes and their separation distance. Thus, spatial variability
in (`\theta`{=tex}), (`\psi`{=tex}), and (K) develops dynamically even
though the underlying soil hydraulic parameterization is uniform with
depth.

DSBM uses adaptive sub-timesteps to limit the amount of water
transferred during a calculation and maintain a stable, mass-conserving
solution.

An optional lookup table can be used for the Clapp-Hornberger hydraulic
functions.

## Noah-MP comparison solver

The repository includes `src/noahmp_soilwater_stateless.c`, an
independent C implementation of the Noah-MP vertical soil-water solution
used for comparison with DSBM.

For the comparison experiments in this repository, Noah-MP is configured
with the same uniform soil hydraulic characteristics throughout the soil
column and the same vertical discretization used by DSBM. The objective
is to compare the numerical soil-water formulations rather than to
compare different soil profiles.

The Noah-MP implementation retains the essential numerical structure of
the original solver:

-   calculation of soil hydraulic diffusivity and conductivity;
-   assembly of the implicit soil-water equations;
-   tridiagonal matrix solution;
-   update of disc soil moisture; and
-   conductivity-based lower-boundary drainage.

The two formulations can therefore be run with the same forcing, soil
hydraulic parameters, discretization, and initial soil-water storage.

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
by the numerical treatment of vertical soil-water movement from
differences caused by the lower-boundary condition.

By default, both comparison implementations use conductivity-based
bottom drainage without imposing field capacity as a percolation
threshold.

The optional switch

``` text
--apply-fc-perc-threshold
```

changes the lower-boundary treatment so that percolation occurs only
from water above field capacity and removal of percolating water does
not reduce the bottom-disc moisture below `theta_fc`.

For Noah-MP, this switch is an experimental comparison option. It is
**not** the native free-drainage behavior represented by the Noah-MP
`OPT_RUN = 3` and `OPT_RUN = 7` soil-water paths associated with the
Schaake and Xinanjiang surface-water partitioning configurations used by
the National Water Model.

The option is retained because it permits controlled comparisons of DSBM
and Noah-MP under either lower-boundary assumption.

## Initial DSBM--Noah-MP comparison

For the current test case, when the field-capacity percolation threshold
is **not** applied, DSBM and the Noah-MP soil-water solution agree very
closely. NSE and KGE values for the simulated soil-moisture states and
water fluxes are generally greater than 0.98.

This result is notable because the two formulations solve the vertical
soil-water problem using substantially different numerical methods:

-   DSBM explicitly transfers water between finite soil discs using
    Darcy-Buckingham fluxes and adaptive sub-timesteps.
-   Noah-MP solves its discretized soil-water equations using an
    implicit tridiagonal matrix solution.

Both calculations in this experiment represent the same vertically
uniform soil hydraulic properties and use comparable boundary
conditions.

The reported NSE and KGE values describe the supplied comparison
experiment. They should not be interpreted as validation across
arbitrary soil properties, forcings, initial conditions, climates, or
numerical configurations.

## Repository layout

``` text
.
├── Makefile
├── configs/
├── include/
│   ├── noahmp_soilwater_stateless.h
│   ├── soil_cli.h
│   ├── soil_config.h
│   ├── soil_data_types.h
│   ├── soil_helpers.h
│   ├── soil_kernel_stateless.h
│   └── util.h
├── src/
│   ├── bmi_soil_driver.c
│   ├── noahmp_soilwater_stateless.c
│   ├── soil_helpers.c
│   ├── soil_kernel_stateless.c
│   └── util.c
├── run_dsbm.csh
├── run_dsbm_fc_limit.csh
├── run_noahmp.csh
├── run_noahmp_fc_limit.csh
├── example_hourly_rainfall_PET.txt
└── 10x_hourly_rainfall_PET.txt
```

Generated executables, build products, and model output are
intentionally excluded from version control.

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

## Running the comparison

Four convenience scripts are supplied:

``` text
run_dsbm.csh
run_noahmp.csh
run_dsbm_fc_limit.csh
run_noahmp_fc_limit.csh
```

The first pair runs the conductivity-based lower-boundary comparison:

``` bash
./run_dsbm.csh
./run_noahmp.csh
```

The second pair applies the optional field-capacity percolation
threshold to both formulations:

``` bash
./run_dsbm_fc_limit.csh
./run_noahmp_fc_limit.csh
```

The scripts use the same uniform soil hydraulic properties, vertical
discretization, forcing, and initial total soil-water storage so that
state and flux time series can be compared directly.

## Driver examples

A DSBM run can be made with:

``` bash
bin/soil_driver \
    --solver dsbm \
    --forcing example_hourly_rainfall_PET.txt \
    --outdir out/dsbm \
    --write-theta \
    --write-fluxes \
    --write-mb \
    --dz 0.1,0.3,0.6,1.0 \
    --zc 0.05,0.25,0.7,1.5 \
    --init-target-storage 0.34
```

For the corresponding Noah-MP calculation:

``` bash
bin/soil_driver \
    --solver noahmp \
    --forcing example_hourly_rainfall_PET.txt \
    --outdir out/noahmp \
    --write-theta \
    --write-fluxes \
    --write-mb \
    --dz 0.1,0.3,0.6,1.0 \
    --zc 0.05,0.25,0.7,1.5 \
    --init-target-storage 0.34
```

DSBM additionally supports a Clapp-Hornberger lookup table:

``` text
--use-lut --lut-n 400 --lut-Theta-min 1e-6
```

## Forcing format

The driver reads hourly forcing records containing:

``` text
YYYY MM DD HH Mi  Precip_mm  PET_mm
```

Lines beginning with `#` are ignored.

Precipitation and potential evapotranspiration are supplied in mm per
hour.

## Output

The driver can write:

-   volumetric soil-moisture time series for each disc;
-   water-flux time series;
-   timestep mass-balance diagnostics; and
-   a mass-balance summary.

Enable these with:

``` text
--write-theta
--write-fluxes
--write-mb
```

The output time coordinate can be an integer timestep index, Julian
Date, or calendar time.

## Soil hydraulic parameters

The driver supports a simple parameter file containing key/value pairs
such as:

``` text
theta_r
theta_sat
theta_fc
theta_wp
theta_aet_eq_pet
K_sat_cm_per_h
phi_sat_cm
b_exp
perc_limiter_0_to_1
klf_m_per_h
```

These are properties or parameters of the **uniform soil column**. The
current DSBM formulation does not assign an independent set of hydraulic
properties to each disc.

## Numerical comparison objective

The DSBM--Noah-MP experiment is designed primarily to answer:

> How closely can a simple discrete soil-moisture balance formulation
> reproduce the vertical soil-moisture states and water fluxes produced
> by the Noah-MP implicit soil-water solver when both represent the same
> uniform soil?

This requires keeping the physical problem as nearly identical as
possible between the two formulations:

-   identical uniform soil hydraulic parameters;
-   identical soil-column depth;
-   identical vertical discretization;
-   identical initial soil-water storage and moisture profile;
-   identical water forcing; and
-   comparable upper and lower boundary conditions.

Differences in the resulting state and flux time series can then be
attributed primarily to the numerical formulations and to any explicitly
selected process differences.

The optional field-capacity percolation threshold provides one
controlled experiment for separating a lower-boundary assumption from
the behavior of the vertical soil-water solution itself.

## Mass conservation

Mass conservation is a primary diagnostic.

Each timestep accounts for precipitation entering the soil, surface
excess, evapotranspiration, lateral flow, bottom percolation, and change
in total soil-water storage.

Comparison experiments should examine the mass-balance residual in
addition to statistical agreement between DSBM and Noah-MP state and
flux time series.

## Status

DSBM is research and development software.

The repository currently contains:

-   the stateless DSBM soil-water formulation;
-   an independent C implementation of the Noah-MP soil-water solution
    for numerical comparison;
-   a common command-line driver;
-   common forcing and configuration capability;
-   mass-balance diagnostics; and
-   reproducible DSBM/Noah-MP comparison scripts.

The initial comparison demonstrates close agreement for the supplied
experiment. Broader evaluation should include different uniform soil
hydraulic properties, initial moisture states, precipitation regimes,
evapotranspiration demand, and numerical timestep conditions.

## Author

Fred L. Ogden\
NOAA/National Weather Service\
Office of Water Prediction

DSBM originated during development of the CFE soil-water formulation and
subsequent investigation of a computationally simple, discretized
representation of vertical soil-moisture dynamics.

## License

Source files originating from CFE/DSBM development contain Apache
License 2.0 notices.

Before assigning a repository-wide license, the provenance and
applicable licensing terms of the Noah-MP comparison implementation
should be documented explicitly.
