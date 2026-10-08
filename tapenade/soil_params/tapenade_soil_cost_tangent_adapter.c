/*
 * tapenade_soil_cost_tangent_adapter.c
 *
 * The only file that sees Tapenade's generated TANGENT of
 * dsbm_soil_cost_from_params() (E10).  It #includes the generated source,
 * so the call is compiled against the real generated definition and any
 * signature mismatch is a compile error.
 *
 * Expected generated signature (Tapenade tangent mode, suffix "d"):
 *
 *   void dsbm_soil_cost_from_params_d(
 *       double klf_m_per_h, double klf_m_per_hd,
 *       double perc_limiter_0_to_1, double perc_limiter_0_to_1d,
 *       double K_sat_cm_per_h, double K_sat_cm_per_hd,
 *       double b_exp, double b_expd,
 *       const double *theta_in, int n_steps,
 *       int n_sub_setting, [int n_sub_settingd,]
 *       const double *rain_mm_per_h, const double *pet_mm_per_h,
 *       const double *observed_lateral_m, const double *observed_percolation_m,
 *       double sst_lateral_m2, double sst_percolation_m2,
 *       double *cost_function_value, double *cost_function_valued);
 *
 * run_tapenade_soil_params.csh defines SOIL_COST_TANGENT_HAS_NSUB_SLOT if
 * Tapenade emits the integer slot n_sub_settingd (integers have no
 * derivative; 0 is passed).
 *
 * One tangent call gives the directional derivative of the cost function
 * along the seed vector (klf, perc_limiter, Ksat, b), in physical units.
 *
 * ASCII only.
 */

#include "tapenade_soil_cost_primal_d.c"

void soil_cost_tangent_from_tapenade(const double parameters[4],
                                     const double seed[4],
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
#ifdef SOIL_COST_TANGENT_HAS_NSUB_SLOT
    const int n_sub_setting_tangent = 0;   /* integer: no derivative */
#endif

    dsbm_soil_cost_from_params_d(parameters[0], seed[0],
                                 parameters[1], seed[1],
                                 parameters[2], seed[2],
                                 parameters[3], seed[3],
                                 (double *)theta_in, n_steps,
                                 n_sub_setting,
#ifdef SOIL_COST_TANGENT_HAS_NSUB_SLOT
                                 n_sub_setting_tangent,
#endif
                                 (double *)rain_mm_per_h, (double *)pet_mm_per_h,
                                 (double *)observed_lateral_m,
                                 (double *)observed_percolation_m,
                                 sst_lateral_m2, sst_percolation_m2,
                                 cost_function_value, directional_derivative_of_cost);
}
