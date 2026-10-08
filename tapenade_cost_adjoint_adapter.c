/*
 * tapenade_cost_adjoint_adapter.c
 *
 * The only file that sees Tapenade's generated REVERSE-mode (adjoint)
 * code for dsbm_cost_from_params().  It #includes the generated source,
 * so the call is compiled against the real generated definition and any
 * signature mismatch is a compile error.
 *
 * Expected generated signature (Tapenade reverse mode, suffix "b"):
 *
 *   void dsbm_cost_from_params_b(
 *       double klf_m_per_h, double *klf_m_per_hb,
 *       double perc_limiter_0_to_1, double *perc_limiter_0_to_1b,
 *       const double *theta_in, int n_steps, int n_sub_setting,
 *       const double *rain_mm_per_h, const double *pet_mm_per_h,
 *       const double *observed_lateral_m, const double *observed_percolation_m,
 *       double sst_lateral_m2, double sst_percolation_m2,
 *       double *cost_function_value, double *cost_function_valueb);
 *
 * Semantics of reverse mode: on entry *cost_function_valueb holds the
 * adjoint seed (1 for a gradient); on exit the parameter adjoints
 * *klf_m_per_hb and *perc_limiter_0_to_1b have been INCREMENTED by
 * seed * d(cost)/d(parameter), so they are zeroed first.  One call gives
 * the whole gradient, whatever the number of parameters.
 *
 * If Tapenade adds an integer adjoint slot for n_sub_setting,
 * run_tapenade_adjoint.csh defines COST_ADJOINT_NSUB_SLOT_BY_VALUE
 * (int n_sub_settingb) or COST_ADJOINT_NSUB_SLOT_BY_POINTER
 * (int *n_sub_settingb).  Integers have no derivative.
 *
 * The generated code uses the ADFirstAidKit stack (adStack.c), which
 * stores the values the reverse sweep needs from the forward sweep.
 *
 * ASCII only.
 */

#include "tapenade_cost_primal_b.c"

void cost_adjoint_from_tapenade(double klf_m_per_h,
                                double perc_limiter_0_to_1,
                                const double *theta_in,
                                int n_steps,
                                int n_sub_setting,
                                const double *rain_mm_per_h,
                                const double *pet_mm_per_h,
                                const double *observed_lateral_m,
                                const double *observed_percolation_m,
                                double sst_lateral_m2,
                                double sst_percolation_m2,
                                double cost_adjoint_seed,
                                double *d_cost_d_klf,
                                double *d_cost_d_perc_limiter)
{
    double cost_value = 0.0;
    double cost_adjoint = cost_adjoint_seed;
#ifdef COST_ADJOINT_NSUB_SLOT_BY_POINTER
    int n_sub_setting_adjoint = 0;
#endif

    *d_cost_d_klf = 0.0;
    *d_cost_d_perc_limiter = 0.0;

    dsbm_cost_from_params_b(klf_m_per_h, d_cost_d_klf,
                            perc_limiter_0_to_1, d_cost_d_perc_limiter,
                            (double *)theta_in, n_steps,
                            n_sub_setting,
#ifdef COST_ADJOINT_NSUB_SLOT_BY_VALUE
                            0,
#endif
#ifdef COST_ADJOINT_NSUB_SLOT_BY_POINTER
                            &n_sub_setting_adjoint,
#endif
                            (double *)rain_mm_per_h, (double *)pet_mm_per_h,
                            (double *)observed_lateral_m,
                            (double *)observed_percolation_m,
                            sst_lateral_m2, sst_percolation_m2,
                            &cost_value, &cost_adjoint);
}
