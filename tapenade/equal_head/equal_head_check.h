/*
 * equal_head_check.h
 *
 * Shared index definitions for the E11 equal-head limiter derivative
 * check (FLO 11/06/2026).  Deliberately includes no generated Tapenade
 * header, so the driver can include the production headers.
 *
 * ASCII only.
 */

#ifndef EQUAL_HEAD_CHECK_H
#define EQUAL_HEAD_CHECK_H

/* all arguments of limit_downward_transfer_to_equal_head(), in order */
#define INPUT_TRANSFER_REQUESTED_M   0
#define INPUT_THETA_UPPER            1
#define INPUT_THETA_LOWER            2
#define INPUT_DZ_UPPER_M             3
#define INPUT_DZ_LOWER_M             4
#define INPUT_THETA_SAT              5
#define INPUT_PHI_SAT_CM             6
#define INPUT_B_EXP                  7
#define N_LIMITER_INPUTS             8

/* the active (differentiated) inputs, in order */
#define ACTIVE_TRANSFER_REQUESTED_M  0
#define ACTIVE_THETA_UPPER           1
#define ACTIVE_THETA_LOWER           2
#define ACTIVE_THETA_SAT             3
#define ACTIVE_PHI_SAT_CM            4
#define ACTIVE_B_EXP                 5
#define N_LIMITER_ACTIVE             6

/* in equal_head_tangent_adapter.c */
void equal_head_tangent_from_tapenade(const double limiter_inputs[N_LIMITER_INPUTS],
                                      const double seed[N_LIMITER_ACTIVE],
                                      double *limited_transfer_m,
                                      double *directional_derivative_m);

#endif /* EQUAL_HEAD_CHECK_H */
