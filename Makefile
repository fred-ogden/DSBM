# =========================
# Soil stateless kernel — Makefile
# =========================

CC      ?= cc
AR      ?= ar
RM      ?= rm -f

# ---- Configuration ----
CSTD        ?= c11
NDISC       ?= 4
THETA_MIN   ?= 1.0e-03
MODE        ?= release

EXTRA_CFLAGS  ?=
EXTRA_LDFLAGS ?=

# ---- Directories ----
INCDIR   := include
SRCDIR   := src
BUILDDIR := build
OBJDIR   := $(BUILDDIR)/obj
LIBDIR   := $(BUILDDIR)/lib
BINDIR   := bin

# ---- Outputs ----
LIBSOIL  := $(LIBDIR)/libsoil.a
DRIVER   := $(BINDIR)/soil_driver

# ---- Sources ----
LIB_SRC  := \
  $(SRCDIR)/soil_helpers.c \
  $(SRCDIR)/dsbm_soilmoisture_stateless.c \
  $(SRCDIR)/noahmp_soilmoisture_stateless.c \
  $(SRCDIR)/util.c

DRV_SRC  := \
  $(SRCDIR)/bmi_soil_driver.c

LIB_OBJ  := $(LIB_SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
DRV_OBJ  := $(DRV_SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
DEPS     := $(LIB_OBJ:.o=.d) $(DRV_OBJ:.o=.d)

WARNFLAGS := -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes
CPPFLAGS  := -I$(INCDIR) -DNDISC=$(NDISC) -DTHETA_MIN=$(THETA_MIN)

CFLAGS_common := -std=$(CSTD) $(WARNFLAGS) -MMD -MP $(CPPFLAGS) $(EXTRA_CFLAGS)

ifeq ($(MODE),debug)
  CFLAGS  := $(CFLAGS_common) -O0 -g
  LDFLAGS := $(EXTRA_LDFLAGS)
else
  CFLAGS  := $(CFLAGS_common) -O3 -DNDEBUG
  LDFLAGS := $(EXTRA_LDFLAGS)
endif

LDLIBS := -lm

# ---- Default ----
.PHONY: all
all: $(DRIVER)

# ---- Link driver ----
$(DRIVER): $(LIBSOIL) $(DRV_OBJ) | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(DRV_OBJ) $(LIBSOIL) $(LDFLAGS) $(LDLIBS)

# ---- Static library ----
$(LIBSOIL): $(LIB_OBJ) | $(LIBDIR)
	$(AR) rcs $@ $(LIB_OBJ)

# ---- Compile ----
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# ---- Dirs ----
$(OBJDIR) $(LIBDIR) $(BINDIR):
	@mkdir -p $@

# ---- Branch-activity census build ----
# bin/soil_driver_census: the same driver and kernel with the DSBM_CENSUS()
# hooks compiled in (-DDSBM_CENSUS_BUILD).  It prints a census of active
# thresholds, caps and branches at exit.  Objects go to a separate
# directory so the normal build is never mixed with census objects.
CENSUSOBJDIR := $(BUILDDIR)/obj_census
CENSUS_SRC   := $(LIB_SRC) $(DRV_SRC) $(SRCDIR)/dsbm_census.c
CENSUS_OBJ   := $(CENSUS_SRC:$(SRCDIR)/%.c=$(CENSUSOBJDIR)/%.o)
CENSUS_BIN   := $(BINDIR)/soil_driver_census

.PHONY: census
census: $(CENSUS_BIN)

$(CENSUS_BIN): $(CENSUS_OBJ) | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(CENSUS_OBJ) $(LDFLAGS) $(LDLIBS)

$(CENSUSOBJDIR)/%.o: $(SRCDIR)/%.c | $(CENSUSOBJDIR)
	$(CC) $(CFLAGS) -DDSBM_CENSUS_BUILD -c $< -o $@

$(CENSUSOBJDIR):
	@mkdir -p $@

# ---- Convenience ----
.PHONY: release debug clean veryclean print
release: ; @$(MAKE) MODE=release
debug:   ; @$(MAKE) MODE=debug

clean:
	$(RM) $(LIB_OBJ) $(DRV_OBJ) $(DEPS)

veryclean: clean
	$(RM) -r $(BUILDDIR) $(BINDIR)

print:
	@echo "CC=$(CC)"
	@echo "MODE=$(MODE)"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "NDISC=$(NDISC)"
	@echo "THETA_MIN=$(THETA_MIN)"
	@echo "LIBSOIL=$(LIBSOIL)"
	@echo "DRIVER=$(DRIVER)"

# ---- Run helper ----
# Usage examples:
#   make run FORCE=example_hourly_rainfall_PET.txt OUT=out
#   make run FORCE=forcing.txt OUT=run TIMECOL=calendar VERB=2
#
# You can override:
#   S0 (target storage, m), DZ, ZC, TIMECOL, VERB, TSTART
#
FORCE    ?=
OUT      ?= out
S0       ?= 0.34
# geometry (comma-separated, as required by the driver parser)
DZ       ?= 0.1,0.3,0.6,1.0
ZC       ?= 0.05,0.25,0.7,1.5
TIMECOL  ?= index
VERB     ?= 1
TSTART   ?= 0

.PHONY: run
run: $(DRIVER)
	@[ -n "$(FORCE)" ] || (echo "Set FORCE=path/to/forcing.txt"; exit 2)
	@mkdir -p $(OUT)
	@$(DRIVER) \
	  --forcing $(FORCE) \
	  --outdir $(OUT) \
	  --verbosity $(VERB) \
	  --write-theta --write-fluxes --write-mb \
	  --timecol $(TIMECOL) --tindex-start $(TSTART) \
	  --use-lut --lut-n 400 --lut-Theta-min 1e-6 \
	  --dz $(DZ) --zc $(ZC) \
	  --init-target-storage $(S0)

# ---- Dependencies ----
-include $(DEPS)
-include $(CENSUS_OBJ:.o=.d)

