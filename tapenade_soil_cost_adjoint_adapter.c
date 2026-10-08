/*
 * tapenade_soil_cost_adjoint_adapter.c
 *
 * The only file that sees Tapenade's generated REVERSE-mode (adjoint) code
 * for dsbm_soil_cost_from_params() (E10).  It #includes the generated
 * source, so the call is compiled against the real generated definition
 * and any signature mismatch is a compile error.
 *
 * Expected generated signature (Tapenade reverse mode, suffix "b"):
 *
 *   void dsbm_soil_cost_from_params_b(
 *       double klf_m_per_h, double *klf_m_per_hb,
 *       double perc_limiter_0_to_1, double *perc_limiter_0_to_1b,
 *       double K_sat_cm_per_h, double *K_sat_cm_per_hb,
 *       double b_exp, double *b_expb,
 *       const double *theta_in, int n_steps,
 *       int n_sub_setting, [int *n_sub_settingb or int n_sub_settingb,]
 *       const double *rain_mm_per_h, const double *pet_mm_per_h,
 *       const double *observed_lateral_m, const double *observed_percolation_m,
 *       double sst_lateral_m2, double sst_percolation_m2,
 *       double *cost_function_value, double *cost_function_valueb);
 *
 * On entry *cost_function_valueb holds the adjoint seed (1 for a
 * gradient); on exit the parameter adjoints have been INCREMENTED by
 * seed * d(cost)/d(parameter), so they are zeroed first.  One call gives
 * all four gradient components, in physical units.
 *
 * run_tapenade_soil_params.csh defines SOIL_COST_ADJOINT_NSUB_SLOT_BY_POINTER
 * or SOIL_COST_ADJOINT_NSUB_SLOT_BY_VALUE if Tapenade emits an integer
 * adjoint slot for n_sub_setting (E9 produced the pointer form).
 *
 * ASCII only.
 */

#include "tapenade_soil_cost_primal_b.c"

void soil_cost_adjoint_from_tapenade(const double parameters[4],
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
                                     double gradient_of_cost[4])
{
    double cost_value = 0.0;
    double cost_adjoint = cost_adjoint_seed;
    int i_parameter;
#ifdef SOIL_COST_ADJOINT_NSUB_SLOT_BY_POINTER
    int n_sub_setting_adjoint = 0;
#endif

    for (i_parameter = 0; i_parameter < 4; i_parameter++) {
        gradient_of_cost[i_parameter] = 0.0;
    }

    dsbm_soil_cost_from_params_b(parameters[0], &gradient_of_cost[0],
                                 parameters[1], &gradient_of_cost[1],
                                 parameters[2], &gradient_of_cost[2],
                                 parameters[3], &gradient_of_cost[3],
                                 (double *)theta_in, n_steps,
                                 n_sub_setting,
#ifdef SOIL_COST_ADJOINT_NSUB_SLOT_BY_VALUE
                                 0,
#endif
#ifdef SOIL_COST_ADJOINT_NSUB_SLOT_BY_POINTER
                                 &n_sub_setting_adjoint,
#endif
                                 (double *)rain_mm_per_h, (double *)pet_mm_per_h,
                                 (double *)observed_lateral_m,
                                 (double *)observed_percolation_m,
                                 sst_lateral_m2, sst_percolation_m2,
                                 &cost_value, &cost_adjoint);
}
