/*
 * dsbm_census.c
 *
 * Branch-activity census for the DSBM kernel.  Compiled only into the
 * census build (make census -> bin/soil_driver_census), where the
 * DSBM_CENSUS() hooks in the kernel call dsbm_census_count().  At program
 * exit it prints how often each threshold, cap, and branch was active.
 *
 * Counts are reported per occurrence, as a percentage of the opportunities
 * for that event (substeps, interface-substeps, disc-substeps, or
 * timesteps), and as the number of timesteps (hours) in which the event
 * occurred at least once.
 *
 * Discrete Soil Moisture Balance Model (DSBM)
 * Author: Fred L. Ogden, Ph.D., P.E., NOAA/National Weather Service
 * SPDX-License-Identifier: Apache-2.0
 *
 * ASCII only.
 */

#include <stdio.h>
#include <stdlib.h>

#include "soil_config.h"
#include "dsbm_census_hooks.h"

#define CENSUS_MAX_INDEX 64   /* indices: discs/interfaces, or n_sub for CENSUS_TIMESTEP */

/* what each event is counted against, for the percentage column */
#define PER_TIMESTEP        0
#define PER_SUBSTEP         1
#define PER_INTERFACE_SUB   2
#define PER_DISC_SUB        3
#define PER_DISC_TIMESTEP   4
#define PER_ROOT_DISC_SUB   5

typedef struct {
    const char *name;
    const char *description;
    int denominator_kind;
    int has_index;            /* 1: report by disc or interface */
} CensusEventInfo;

static const CensusEventInfo census_event_info[CENSUS_N_EVENTS] = {
    {"timestep",          "hourly timesteps (index = n_sub used)",          PER_TIMESTEP,      0},
    {"fc_cross_down",     "disc crossed theta_fc downward over the hour",   PER_DISC_TIMESTEP, 1},
    {"fc_cross_up",       "disc crossed theta_fc upward over the hour",     PER_DISC_TIMESTEP, 1},
    {"wp_cross_down",     "disc crossed theta_wp downward over the hour",   PER_DISC_TIMESTEP, 1},
    {"wp_cross_up",       "disc crossed theta_wp upward over the hour",     PER_DISC_TIMESTEP, 1},
    {"substep",           "substeps",                                       PER_SUBSTEP,       0},
    {"rain",              "substeps with rain",                             PER_SUBSTEP,       0},
    {"rain_capped",       "rain exceeded disc-1 space (surface excess)",    PER_SUBSTEP,       0},
    {"flux_down",         "Darcy-Buckingham flux downward",                 PER_INTERFACE_SUB, 1},
    {"flux_up",           "Darcy-Buckingham flux upward",                   PER_INTERFACE_SUB, 1},
    {"chain_cap",         "downward pass limited by downstream acceptance", PER_INTERFACE_SUB, 1},
    {"down_cap_max_out",  "downward transfer limited by max_out",           PER_INTERFACE_SUB, 1},
    {"down_cap_donor",    "downward transfer limited by donor floor",       PER_INTERFACE_SUB, 1},
    {"down_cap_receiver", "downward transfer limited by receiver space",    PER_INTERFACE_SUB, 1},
    {"up_cap_donor",      "upward transfer limited by donor floor",         PER_INTERFACE_SUB, 1},
    {"up_cap_receiver",   "upward transfer limited by receiver space",      PER_INTERFACE_SUB, 1},
    {"transfer_clamp",    "post-transfer clamp changed theta (roundoff)",   PER_INTERFACE_SUB, 1},
    {"perc_active",       "bottom potential percolation > 0",               PER_SUBSTEP,       0},
    {"perc_capped",       "percolation limited by disc-4 water",            PER_SUBSTEP,       0},
    {"et_dry",            "ET disc at or below theta_wp (no AET)",          PER_ROOT_DISC_SUB, 1},
    {"et_stressed",       "ET disc between theta_wp and theta_fc",          PER_ROOT_DISC_SUB, 1},
    {"et_unstressed",     "ET disc at or above theta_fc (AET = PET)",       PER_ROOT_DISC_SUB, 1},
    {"et_capped",         "AET limited by water above theta_wp",            PER_ROOT_DISC_SUB, 1},
    {"lat_active",        "disc above theta_fc at lateral step",            PER_DISC_SUB,      1},
    {"lat_euler_cap",     "forward-Euler lateral limited (legacy only)",    PER_DISC_SUB,      1},
    {"safety_clamp",      "end-of-substep safety clamp changed theta",      PER_DISC_SUB,      1}
};

static long long census_count[CENSUS_N_EVENTS][CENSUS_MAX_INDEX];
static long long census_hours[CENSUS_N_EVENTS][CENSUS_MAX_INDEX];
static long long census_last_hour[CENSUS_N_EVENTS][CENSUS_MAX_INDEX];
static long long census_current_hour = 0;
static long long census_total_substeps = 0;
static int census_registered = 0;

static void dsbm_census_report(void);

void dsbm_census_count(int event_id, int index)
{
    if (!census_registered) {
        for (int e = 0; e < CENSUS_N_EVENTS; e++) {
            for (int k = 0; k < CENSUS_MAX_INDEX; k++) {
                census_last_hour[e][k] = -1;
            }
        }
        atexit(dsbm_census_report);
        census_registered = 1;
    }
    if (event_id < 0 || event_id >= CENSUS_N_EVENTS) return;
    if (index < 0) index = 0;
    if (index >= CENSUS_MAX_INDEX) index = CENSUS_MAX_INDEX - 1;

    if (event_id == CENSUS_TIMESTEP) {
        census_current_hour++;
        census_total_substeps = census_total_substeps + index;
    }

    census_count[event_id][index]++;
    if (census_last_hour[event_id][index] != census_current_hour) {
        census_hours[event_id][index]++;
        census_last_hour[event_id][index] = census_current_hour;
    }
}

static double census_denominator(int kind)
{
    double n_hours = (double)census_current_hour;
    double n_sub = (double)census_total_substeps;
    double denominator;

    if (kind == PER_TIMESTEP) {
        denominator = n_hours;
    } else if (kind == PER_SUBSTEP) {
        denominator = n_sub;
    } else if (kind == PER_INTERFACE_SUB) {
        denominator = n_sub;        /* per interface: each interface once per substep */
    } else if (kind == PER_DISC_SUB || kind == PER_ROOT_DISC_SUB) {
        denominator = n_sub;        /* per disc: each disc once per substep */
    } else {
        denominator = n_hours;      /* per disc, per timestep */
    }
    if (denominator < 1.0) denominator = 1.0;
    return denominator;
}

static void dsbm_census_report(void)
{
    double n_hours = (double)census_current_hour;
    int n_index_shown;

    printf("\n================ DSBM BRANCH-ACTIVITY CENSUS ================\n");
    if (n_hours < 1.0) n_hours = 1.0;
    printf("timesteps (hours) = %lld   substeps = %lld   mean substeps/hour = %.3f\n",
           census_current_hour, census_total_substeps,
           (double)census_total_substeps / n_hours);

    printf("\nn_sub used per hour:\n");
    for (int k = 1; k < CENSUS_MAX_INDEX; k++) {
        if (census_count[CENSUS_TIMESTEP][k] > 0) {
            printf("   n_sub = %2d : %10lld hours (%7.3f %%)\n", k,
                   census_count[CENSUS_TIMESTEP][k],
                   100.0 * (double)census_count[CENSUS_TIMESTEP][k] / n_hours);
        }
    }

    printf("\nPercent of opportunities in which the event occurred, and the\n");
    printf("number of hours in which it occurred at least once.  Interface i\n");
    printf("is between disc i and disc i+1 (discs numbered 1..%d from the top).\n", NDISC);

    for (int e = 1; e < CENSUS_N_EVENTS; e++) {
        const CensusEventInfo *info = &census_event_info[e];
        double denominator = census_denominator(info->denominator_kind);

        if (info->has_index) {
            n_index_shown = NDISC;
            if (info->denominator_kind == PER_INTERFACE_SUB) n_index_shown = NDISC - 1;
            printf("\n%-18s %s\n", info->name, info->description);
            for (int k = 0; k < n_index_shown; k++) {
                const char *label = "disc";
                if (info->denominator_kind == PER_INTERFACE_SUB) label = "intf";
                printf("   %s %d : %12lld  %9.4f %%  in %8lld hours\n",
                       label, k + 1, census_count[e][k],
                       100.0 * (double)census_count[e][k] / denominator,
                       census_hours[e][k]);
            }
        } else {
            printf("\n%-18s %s\n", info->name, info->description);
            printf("            %12lld  %9.4f %%  in %8lld hours\n",
                   census_count[e][0],
                   100.0 * (double)census_count[e][0] / denominator,
                   census_hours[e][0]);
        }
    }
    printf("=============================================================\n");
}
