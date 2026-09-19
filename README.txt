Noah-MP / DSBM comparison patch

Files:
  src/noahmp_soilwater_stateless.c
  include/noahmp_soilwater_stateless.h
  src/bmi_soil_driver.c          modified: --solver dsbm|noahmp
  include/soil_cli.h             modified: stores solver selection

Build must include src/noahmp_soilwater_stateless.c in addition to the
existing DSBM sources.  NDISC remains a compile-time definition as in the
existing project.

Example:
  ./soil_driver --solver noahmp ...
  ./soil_driver --solver dsbm ...

The Noah-MP comparison kernel translates the supplied WDFCND1/SRT/SSTEP/
ROSR12 path for unfrozen soil.  It intentionally does not import the full
Noah-MP runoff/groundwater machinery.

Comparison-wrapper choices are documented at the top of
src/noahmp_soilwater_stateless.c.  In particular, DWSAT is derived from the
same Clapp-Hornberger relations used by DSBM because the supplied Noah-MP
files use parameters%DWSAT but do not contain its parameter-construction
code.
