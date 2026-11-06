/*
 * equal_head_tangent_adapter.c
 *
 * The only file that sees Tapenade's generated TANGENT of
 * equal_head_check_primal() (E11, FLO 11/06/2026).  It #includes the
 * generated source, so the call is compiled against the real generated
 * definition and any signature mismatch is a compile error.
 *
 * Expected generated signature (Tapenade tangent mode, suffix "d"):
 *
 *   void equal_head_check_primal_d(
 *       double transfer_requested_m, double transfer_requested_md,
 *       double theta_upper, double theta_upperd,
 *       double theta_lower, double theta_lowerd,
 *       double dz_upper_m, double dz_lower_m,
 *       double theta_sat, double theta_satd,
 *       double phi_sat_cm, double phi_sat_cmd,
 *       double b_exp, double b_expd,
 *       double *limited_transfer_m, double *limited_transfer_md);
 *
 * One call gives the directional derivative of the limited transfer
 * along the seed vector over the six active inputs, in physical units.
 *
 * ASCII only.
 */

#include "equal_head_check_primal_d.c"

#include "equal_head_check.h"

void equal_head_tangent_from_tapenade(const double limiter_inputs[N_LIMITER_INPUTS],
                                      const double seed[N_LIMITER_ACTIVE],
                                      double *limited_transfer_m,
                                      double *directional_derivative_m)
{
    equal_head_check_primal_d(limiter_inputs[INPUT_TRANSFER_REQUESTED_M], seed[ACTIVE_TRANSFER_REQUESTED_M],
                              limiter_inputs[INPUT_THETA_UPPER],          seed[ACTIVE_THETA_UPPER],
                              limiter_inputs[INPUT_THETA_LOWER],          seed[ACTIVE_THETA_LOWER],
                              limiter_inputs[INPUT_DZ_UPPER_M],
                              limiter_inputs[INPUT_DZ_LOWER_M],
                              limiter_inputs[INPUT_THETA_SAT],            seed[ACTIVE_THETA_SAT],
                              limiter_inputs[INPUT_PHI_SAT_CM],           seed[ACTIVE_PHI_SAT_CM],
                              limiter_inputs[INPUT_B_EXP],                seed[ACTIVE_B_EXP],
                              limited_transfer_m, directional_derivative_m);
}
