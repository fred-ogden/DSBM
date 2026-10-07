/*
 * dsbm_census_hooks.h
 *
 * Branch-activity census hooks for the DSBM kernel.
 *
 * In the normal build (and in Tapenade's preprocessing) DSBM_CENSUS(...)
 * expands to nothing, so the kernel is unchanged.  Only the census build
 * (make census, which defines DSBM_CENSUS_BUILD) turns each hook into a
 * call to dsbm_census_count(), which counts how often each threshold,
 * cap, or branch is active.  See src/dsbm_census.c and
 * run_branch_census.csh.
 *
 * Event identifiers are plain macros so they vanish with the hooks.
 * The second hook argument is a disc or interface index (0-based), or 0.
 *
 * Discrete Soil Moisture Balance Model (DSBM)
 * Author: Fred L. Ogden, Ph.D., P.E., NOAA/National Weather Service
 * SPDX-License-Identifier: Apache-2.0
 *
 * ASCII only.
 */

#ifndef DSBM_CENSUS_HOOKS_H
#define DSBM_CENSUS_HOOKS_H

/* ---- once per timestep ---- */
#define CENSUS_TIMESTEP              0   /* index = n_sub used this timestep */
#define CENSUS_FC_CROSS_DOWN         1   /* disc went from > theta_fc to <= theta_fc */
#define CENSUS_FC_CROSS_UP           2   /* disc went from <= theta_fc to > theta_fc */
#define CENSUS_WP_CROSS_DOWN         3   /* disc went from > theta_wp to <= theta_wp */
#define CENSUS_WP_CROSS_UP           4   /* disc went from <= theta_wp to > theta_wp */

/* ---- once per substep ---- */
#define CENSUS_SUBSTEP               5
#define CENSUS_RAIN                  6   /* rain this substep */
#define CENSUS_RAIN_CAPPED           7   /* rain exceeded disc-0 storage space: surface excess */

/* ---- per interface, per substep (index = upper disc of the interface) ---- */
#define CENSUS_FLUX_DOWN             8   /* Darcy-Buckingham flux downward */
#define CENSUS_FLUX_UP               9   /* Darcy-Buckingham flux upward */
#define CENSUS_CHAIN_CAP            10   /* downward pass limited by downstream acceptance */
#define CENSUS_DOWN_CAP_MAX_OUT     11   /* downward transfer limited by max_out */
#define CENSUS_DOWN_CAP_DONOR       12   /* downward transfer limited by donor water above floor */
#define CENSUS_DOWN_CAP_RECEIVER    13   /* downward transfer limited by receiver space to saturation */
#define CENSUS_UP_CAP_DONOR         14   /* upward transfer limited by donor water above floor */
#define CENSUS_UP_CAP_RECEIVER      15   /* upward transfer limited by receiver space to saturation */
#define CENSUS_TRANSFER_CLAMP       16   /* post-transfer floor/saturation clamp changed theta */

/* ---- bottom boundary, per substep ---- */
#define CENSUS_PERC_ACTIVE          17   /* bottom potential percolation > 0 */
#define CENSUS_PERC_CAPPED          18   /* percolation limited by water available in disc 4 */

/* ---- ET, per root disc, per substep (index = disc) ---- */
#define CENSUS_ET_DRY               19   /* theta <= theta_wp: no AET */
#define CENSUS_ET_STRESSED          20   /* theta_wp < theta < theta_aet_eq_pet: linear stress */
#define CENSUS_ET_UNSTRESSED        21   /* theta >= theta_aet_eq_pet: AET = assigned PET */
#define CENSUS_ET_CAPPED            22   /* AET limited by water above theta_wp */

/* ---- lateral, per disc, per substep (index = disc) ---- */
#define CENSUS_LAT_ACTIVE           23   /* disc above theta_fc at the lateral step */
#define CENSUS_LAT_EULER_CAP        24   /* forward-Euler removal limited by water above theta_fc */

/* ---- end-of-substep safety clamp, per disc ---- */
#define CENSUS_SAFETY_CLAMP         25   /* theta changed by the floor/saturation safety clamp */

#define CENSUS_N_EVENTS             26

/*
 * The hook supplies its own semicolon and is written in the kernel
 * WITHOUT one, so the normal build is left with no tokens at all: no
 * empty statements and no empty if or loop bodies.  (Empty bodies are
 * the suspected cause of a Tapenade 3.16 crash in flow-graph
 * differentiation, seen after the hooks were added.)  Code that
 * exists only for the census, such as an if whose only content is a
 * hook, is placed inside #ifdef DSBM_CENSUS_BUILD instead.
 */
#ifdef DSBM_CENSUS_BUILD
void dsbm_census_count(int event_id, int index);
#define DSBM_CENSUS(event_id, index) dsbm_census_count((event_id), (index));
#else
#define DSBM_CENSUS(event_id, index)
#endif

#endif /* DSBM_CENSUS_HOOKS_H */
