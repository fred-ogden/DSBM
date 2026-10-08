/*
 * tapenade_klf_tangent_adapter.c
 *
 * The only file that sees Tapenade's generated tangent of
 * dsbm_lateral_from_klf().  It #includes the generated source directly,
 * so the call below is compiled against the REAL generated definition:
 * a wrong argument count or a wrong argument type is a compile error,
 * not a silent pointer shift like the earlier five-argument mistake in
 * the Ksat experiment.
 *
 * Because this file includes Tapenade's generated *_d.h headers (which
 * repeat the primal typedefs), it must NOT include the production
 * soil_data_types.h.  The test driver includes the production headers
 * and calls only klf_tangent_from_tapenade(), whose signature is fixed
 * here and does not depend on Tapenade.
 *
 * Expected generated signature (Tapenade places each tangent argument
 * immediately after its primal argument, suffix "d"):
 *
 *   void dsbm_lateral_from_klf_d(
 *       double klf_m_per_h, double klf_m_per_hd,
 *       const double *theta_in, int n_steps,
 *       int n_sub_fixed, int n_sub_fixedd,
 *       int lateral_analytic,
 *       int n_sub_minimum, int n_sub_minimumd,
 *       int lateral_substep_severity, int lateral_substep_severityd,
 *       const double *rain_mm_per_h, const double *pet_mm_per_h,
 *       double *lateral_total_by_disc_m, double *lateral_total_by_disc_md,
 *       double *theta_out, double *theta_outd,
 *       double *percolation_total_m, double *percolation_total_md,
 *       double *aet_total_m, double *aet_total_md,
 *       double *precipitation_excess_total_m,
 *       double *precipitation_excess_total_md,
 *       int *n_sub_used_by_step);
 *
 * n_sub_fixedd, n_sub_minimumd, lateral_substep_severityd: the wrapper
 * copies these integers into the SoilControl structure, which Tapenade
 * differentiates as a whole, so it emits a tangent slot for each.
 * lateral_analytic is only tested in an if statement in the wrapper, so
 * it gets no tangent slot.  Integer switches have no derivative; the
 * adapter passes 0 for every integer tangent slot.
 *
 * run_tapenade_klf.csh prints the actual generated prototype so it can be
 * compared with this.
 *
 * ASCII only.
 */

#include "tapenade_klf_primal_d.c"

void klf_tangent_from_tapenade(double klf_m_per_h,
                               const double *theta_in,
                               int n_steps,
                               int n_sub_fixed,
                               int lateral_analytic,
                               int n_sub_minimum,
                               int lateral_substep_severity,
                               const double *rain_mm_per_h,
                               const double *pet_mm_per_h,
                               double *lateral_total_by_disc_m,
                               double *d_lateral_total_by_disc_m_d_klf,
                               double *theta_out,
                               double *d_theta_out_d_klf,
                               double *percolation_total_m,
                               double *d_percolation_total_m_d_klf,
                               double *aet_total_m,
                               double *d_aet_total_m_d_klf,
                               double *precipitation_excess_total_m,
                               double *d_precipitation_excess_total_m_d_klf,
                               int *n_sub_used_by_step)
{
    /* seed: d(klf)/d(klf) = 1 */
    const double klf_seed = 1.0;
    /* integer switches: no tangent */
    const int n_sub_fixed_tangent = 0;
    const int n_sub_minimum_tangent = 0;
    const int lateral_substep_severity_tangent = 0;

    dsbm_lateral_from_klf_d(klf_m_per_h, klf_seed,
                            theta_in, n_steps, n_sub_fixed, n_sub_fixed_tangent,
                            lateral_analytic,
                            n_sub_minimum, n_sub_minimum_tangent,
                            lateral_substep_severity, lateral_substep_severity_tangent,
                            rain_mm_per_h, pet_mm_per_h,
                            lateral_total_by_disc_m, d_lateral_total_by_disc_m_d_klf,
                            theta_out, d_theta_out_d_klf,
                            percolation_total_m, d_percolation_total_m_d_klf,
                            aet_total_m, d_aet_total_m_d_klf,
                            precipitation_excess_total_m,
                            d_precipitation_excess_total_m_d_klf,
                            n_sub_used_by_step);
}
