
#include "soil_helpers.h"

// FLO 10/26: Isolate the active equal-head limiter so Tapenade
// differentiates the production Newton solver, rather than an
// inactive limiter branch in a longer DSBM simulation.
void equal_head_ad_primal(double b_exp, double *limited_transfer_m)
{
    *limited_transfer_m = limit_downward_transfer_to_equal_head(
        0.015, 0.380, 0.200, 0.10, 0.30,
        0.439, 7.275, b_exp);
}
