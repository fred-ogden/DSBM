
#include <math.h>
#include <stdio.h>

void equal_head_ad_primal_b(double b_exp, double *b_expb,
                            double *limited_transfer_m,
                            double *limited_transfer_mb);

int main(void)
{
    const double b_exp = 5.5;
    const double expected_derivative = -2.679234718888122e-4;
    double b_expb = 0.0;
    double limited_transfer_m = 0.0;
    double limited_transfer_mb = 1.0;
    double relative_error;

    // FLO 10/26: Verify reverse-mode differentiation through the
    // active equal-head Newton solver against the tangent derivative.
    equal_head_ad_primal_b(
        b_exp, &b_expb, &limited_transfer_m, &limited_transfer_mb);

    relative_error = fabs(b_expb - expected_derivative) /
                     fabs(expected_derivative);

    printf("Tangent derivative: %.15e m\n", expected_derivative);
    printf("Adjoint derivative: %.15e m\n", b_expb);
    printf("Relative difference: %.15e\n", relative_error);

    if (relative_error > 1.0e-10) {
        fprintf(stderr, "FAIL: tangent-adjoint mismatch\n");
        return 1;
    }

    printf("PASS: active-limiter adjoint derivative\n");
    return 0;
}
