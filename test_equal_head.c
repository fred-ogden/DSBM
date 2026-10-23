
#include <math.h>
#include <stdio.h>

#include "soil_helpers.h"

int main(void)
{
    const double theta_upper = 0.380;
    const double theta_lower = 0.200;
    const double dz_upper_m = 0.10;
    const double dz_lower_m = 0.30;
    const double theta_sat = 0.439;
    const double phi_sat_cm = 7.275;
    const double b_exp = 5.5;
    const double requested_m = 0.015;
    const double distance_m = 0.5 * (dz_upper_m + dz_lower_m);

    double limited_m;
    double theta_upper_after;
    double theta_lower_after;
    double psi_upper_m;
    double psi_lower_m;
    double residual_m;

    // FLO 10/26: Exercise the active equal-head limiter directly.
    // Verify that it reduces an overshooting downward transfer and
    // produces equal total hydraulic head between adjacent discs.
    limited_m = limit_downward_transfer_to_equal_head(
        requested_m, theta_upper, theta_lower,
        dz_upper_m, dz_lower_m, theta_sat, phi_sat_cm, b_exp);

    theta_upper_after = theta_upper - limited_m / dz_upper_m;
    theta_lower_after = theta_lower + limited_m / dz_lower_m;

    psi_upper_m = psi_from_theta(
        theta_upper_after, theta_sat, phi_sat_cm, b_exp);
    psi_lower_m = psi_from_theta(
        theta_lower_after, theta_sat, phi_sat_cm, b_exp);

    residual_m = psi_upper_m - psi_lower_m - distance_m;

    printf("Requested transfer: %.15e m\n", requested_m);
    printf("Limited transfer:   %.15e m\n", limited_m);
    printf("Upper theta after:  %.15e\n", theta_upper_after);
    printf("Lower theta after:  %.15e\n", theta_lower_after);
    printf("Equal-head residual: %.15e m\n", residual_m);

    if (!(limited_m > 0.0 && limited_m < requested_m)) {
        fprintf(stderr, "FAIL: limiter did not activate\n");
        return 1;
    }

    if (fabs(residual_m) > 1.0e-10) {
        fprintf(stderr, "FAIL: equal-head residual exceeds tolerance\n");
        return 1;
    }

    printf("PASS: active limiter and equal-head condition\n");
    return 0;
}
