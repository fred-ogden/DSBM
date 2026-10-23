
#include <math.h>
#include <stdio.h>

#include "soil_helpers.h"

void equal_head_ad_primal_d(double b_exp, double b_expd,
                            double *limited_transfer_m,
                            double *limited_transfer_md);

static double primal_transfer(double b_exp)
{
    return limit_downward_transfer_to_equal_head(
        0.015, 0.380, 0.200, 0.10, 0.30,
        0.439, 7.275, b_exp);
}

int main(void)
{
    const double b_exp = 5.5;
    const double relative_steps[] = {
        1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6
    };

    double limited_transfer_m;
    double derivative_ad;
    double derivative_fd;
    double step;
    double relative_error;
    int i;

    // FLO 10/26: Validate the tangent derivative through the active
    // equal-head Newton solver against central finite differences.
    equal_head_ad_primal_d(
        b_exp, 1.0, &limited_transfer_m, &derivative_ad);

    printf("Limited transfer: %.15e m\n", limited_transfer_m);
    printf("AD derivative:    %.15e m\n", derivative_ad);
    printf("\nRelative step    Finite difference       Relative error\n");

    for (i = 0; i < 5; i++) {
        step = relative_steps[i] * b_exp;

        derivative_fd = (
            primal_transfer(b_exp + step) -
            primal_transfer(b_exp - step)
        ) / (2.0 * step);

        relative_error = fabs(derivative_fd - derivative_ad) /
                         fmax(fabs(derivative_ad), 1.0e-15);

        printf("%12.0e    %+.15e    %.3e\n",
               relative_steps[i], derivative_fd, relative_error);
    }

    return 0;
}
