/*
 * equal_head_check_primal.c
 *
 * Primal wrapper for the E11 experiment (FLO 11/06/2026): Tapenade tangent
 * of the equal-head limit on downward transfers, checked against the
 * analytic implicit-function derivative at every limiter activation on
 * the forcing record (see equal_head_check_driver.c).
 *
 * Head function handed to Tapenade:
 *
 *   independents:  transfer_requested_m  theta_upper  theta_lower
 *                  theta_sat  phi_sat_cm  b_exp
 *   dependent:     limited_transfer_m
 *
 * The disc thicknesses are passive (geometry is never calibrated).
 * Ksat enters the limiter only through phi_sat (CFE3.1 regression
 * phi_sat = 10.415 Ksat^-0.3266), so the phi_sat derivative covers Ksat.
 *
 * This wrapper adds no arithmetic: it only passes the result out through
 * a pointer, so that the generated tangent has a predictable signature.
 *
 * ASCII only.
 */

#include "soil_helpers.h"

void equal_head_check_primal(double transfer_requested_m,
                             double theta_upper,
                             double theta_lower,
                             double dz_upper_m,
                             double dz_lower_m,
                             double theta_sat,
                             double phi_sat_cm,
                             double b_exp,
                             double *limited_transfer_m)
{
    *limited_transfer_m = limit_downward_transfer_to_equal_head(
        transfer_requested_m, theta_upper, theta_lower,
        dz_upper_m, dz_lower_m, theta_sat, phi_sat_cm, b_exp);
}
