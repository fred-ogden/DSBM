/*
 * tapenade_cost_tangent_adapter.c
 *
 * The only file that sees Tapenade's generated tangent of
 * dsbm_cost_from_params().  It #includes the generated source, so the
 * call is compiled against the real generated definition and any
 * signature mismatch is a compile error.
 *
 * Expected generated signature:
 *
 *   void dsbm_cost_from_params_d(
 *       double klf_m_per_h, double klf_m_per_hd,
 *       double perc_limiter_0_to_1, double perc_limiter_0_to_1d,
 *       const double *theta_in, int n_steps,
 *       int n_sub_setting, [int n_sub_settingd,]
 *       const double *rain_mm_per_h, const double *pet_mm_per_h,
 *       const double *observed_lateral_m, const double *observed_percolation_m,
 *       double sst_lateral_m2, double sst_percolation_m2,
 *       double *cost_function_value, double *cost_function_valued);
 *
 * Whether Tapenade emits the integer tangent slot n_sub_settingd depends
 * on how it treats the integer copied into SoilControl.  run_tapenade_cost.csh
 * looks for it in the generated prototype and defines
 * COST_TANGENT_HAS_NSUB_SLOT when present.  Integer arguments have no
 * derivative; the adapter passes 0 for that slot.
 *
 * One tangent call gives the directional derivative of the cost function
 * along (seed_klf, seed_perc).  Seeds (1,0) and (0,1) give the two
 * partial derivatives.
 *
 * ASCII only.
 */

#include "tapenade_cost_primal_d.c"

void cost_tangent_from_tapenade(double klf_m_per_h,
                                double perc_limiter_0_to_1,
                                double seed_klf,
                                double seed_perc_limiter,
                                const double *theta_in,
                                int n_steps,
                                int n_sub_setting,
                                const double *rain_mm_per_h,
                                const double *pet_mm_per_h,
                                const double *observed_lateral_m,
                                const double *observed_percolation_m,
                                double sst_lateral_m2,
                                double sst_percolation_m2,
                                double *cost_function_value,
                                double *directional_derivative_of_cost)
{
#ifdef COST_TANGENT_HAS_NSUB_SLOT
    const int n_sub_setting_tangent = 0;   /* integer: no derivative */
#endif

    dsbm_cost_from_params_d(klf_m_per_h, seed_klf,
                            perc_limiter_0_to_1, seed_perc_limiter,
                            (double *)theta_in, n_steps,
                            n_sub_setting,
#ifdef COST_TANGENT_HAS_NSUB_SLOT
                            n_sub_setting_tangent,
#endif
                            (double *)rain_mm_per_h, (double *)pet_mm_per_h,
                            (double *)observed_lateral_m,
                            (double *)observed_percolation_m,
                            sst_lateral_m2, sst_percolation_m2,
                            cost_function_value, directional_derivative_of_cost);
}
