/*
 * tapenade_adjoint_driver.c
 *
 * E9: the adjoint.  Tapenade REVERSE mode on the same cost-function wrapper
 * as E8 (dsbm_cost_from_params in tapenade_cost_primal.c), compared with
 * the tangent (forward) mode validated in E8.
 *
 * Same synthetic-truth problem as E8: observations from the model at
 * k_lf = 2.0e-3 m/h, perc_limiter = 0.5, fixed 4 substeps; cost = sum of
 * (1 - NSE) for hourly lateral flow and percolation.
 *
 *   E9a  Adjoint gradient against tangent gradient at the E8a test point
 *        for fixed 4, 12 and 24 substeps, plus the dot-product test
 *        w * (J v) = (J^T w) . v for several directions v.
 *   E9b  Run time: primal, one tangent direction, one adjoint (whole
 *        gradient), for 1 year of hourly forcing.
 *   E9c  Memory: growth of peak resident memory during the first adjoint
 *        run, for 1 year and for the full 78168-hour record.  The adjoint
 *        stores forward-sweep values on the ADFirstAidKit stack.
 *   E9d  BFGS calibration (E8c start 5e-4 m/h, 0.8) with the adjoint
 *        gradient and with the tangent gradient; same path, compared time.
 *
 * Usage:
 *   tapenade_adjoint_test  output_directory  forcing_csv_file  [experiments]
 *   experiments: letters, default "abcd".
 *
 * Terminology: soil moisture; discs; cost function.
 * ASCII only.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/resource.h>

#include "soil_data_types.h"
#include "soil_helpers.h"


/* ---- functions in tapenade_cost_primal.c and the two adapters ---- */

double cost_experiment_phi_sat_cm(void);

void dsbm_cost_from_params(double klf_m_per_h,
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
                           double *cost_function_value);

void dsbm_series_from_params(double klf_m_per_h,
                             double perc_limiter_0_to_1,
                             const double *theta_in,
                             int n_steps,
                             int n_sub_setting,
                             const double *rain_mm_per_h,
                             const double *pet_mm_per_h,
                             double *lateral_series_m,
                             double *percolation_series_m);

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
                                double *directional_derivative_of_cost);


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
                                double *d_cost_d_perc_limiter);


/* ---- experiment constants ---- */

#define MAX_HOURS                   80000
#define HOURS_ONE_YEAR              8760
#define KLF_TRUE_M_PER_H            2.0e-3
#define PERC_LIMITER_TRUE           0.5
#define N_SUB_FOR_OBSERVATIONS      4
#define PERC_LIMITER_LOWER_BOUND    0.01
#define PERC_LIMITER_UPPER_BOUND    0.999
#define KLF_LOWER_BOUND_M_PER_H     1.0e-5
#define KLF_UPPER_BOUND_M_PER_H     1.0

static const double disc_center_depth_m[NDISC] = {0.05, 0.25, 0.70, 1.50};


/* ---- forcing (whole record) and the current calibration problem ---- */

static double all_rain_mm_per_h[MAX_HOURS];
static double all_pet_mm_per_h[MAX_HOURS];
static int n_hours_available = 0;

typedef struct {
    double theta_in[NDISC];
    int n_steps;
    double *rain_mm_per_h;
    double *pet_mm_per_h;
    double observed_lateral_m[MAX_HOURS];
    double observed_percolation_m[MAX_HOURS];
    double sst_lateral_m2;
    double sst_percolation_m2;
} CalibrationProblem;

static CalibrationProblem problem;

/* peak resident memory before and after the very first adjoint run (1 year) */
static long rss_before_first_adjoint_kb = 0;
static long rss_after_first_adjoint_kb = 0;
static double first_adjoint_seconds = 0.0;

/* 0 = tangent (one run per parameter), 1 = adjoint (one run per gradient) */
static int gradient_method = 0;


static double wall_clock_seconds(void)
{
    struct timespec now;

    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double)now.tv_sec + 1.0e-9 * (double)now.tv_nsec;
}

static long peak_resident_memory_kilobytes(void)
{
    struct rusage usage;

    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;     /* kilobytes on Linux */
}

static const char *gradient_method_name(void)
{
    if (gradient_method == 1) return "adjoint";
    return "tangent";
}

static int read_forcing_csv(const char *path, int n_wanted,
                            double *rain_mm_per_h, double *pet_mm_per_h)
{
    FILE *fp;
    char line[512];
    int n_read = 0;

    fp = fopen(path, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: cannot open forcing file %s\n", path);
        return -1;
    }
    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return -1;
    }
    while (n_read < n_wanted && fgets(line, sizeof(line), fp) != NULL) {
        char *first_comma = strchr(line, ',');
        double rain_value;
        double pet_value;

        if (first_comma == NULL) continue;
        if (sscanf(first_comma + 1, "%lf,%lf", &rain_value, &pet_value) != 2) continue;
        rain_mm_per_h[n_read] = rain_value;
        pet_mm_per_h[n_read] = pet_value;
        n_read++;
    }
    fclose(fp);
    return n_read;
}

static double sum_of_squared_deviations(const double *series, int n)
{
    double mean_value = 0.0;
    double sum_squares = 0.0;
    int i;

    for (i = 0; i < n; i++) mean_value = mean_value + series[i];
    mean_value = mean_value / (double)n;
    for (i = 0; i < n; i++) {
        double deviation = series[i] - mean_value;
        sum_squares = sum_squares + deviation * deviation;
    }
    return sum_squares;
}

/* Set up the synthetic-truth problem over the first n_hours of forcing. */
static void set_up_problem(int n_hours)
{
    double water_table_depth_m;

    problem.n_steps = n_hours;
    problem.rain_mm_per_h = all_rain_mm_per_h;
    problem.pet_mm_per_h = all_pet_mm_per_h;

    initialize_hydrostatic_from_storage(2.0, 0.439, cost_experiment_phi_sat_cm(), 4.05,
                                        disc_center_depth_m, 0.7843498367003,
                                        &water_table_depth_m, problem.theta_in);

    dsbm_series_from_params(KLF_TRUE_M_PER_H, PERC_LIMITER_TRUE, problem.theta_in,
                            problem.n_steps, N_SUB_FOR_OBSERVATIONS,
                            problem.rain_mm_per_h, problem.pet_mm_per_h,
                            problem.observed_lateral_m, problem.observed_percolation_m);

    problem.sst_lateral_m2 = sum_of_squared_deviations(problem.observed_lateral_m, problem.n_steps);
    problem.sst_percolation_m2 = sum_of_squared_deviations(problem.observed_percolation_m, problem.n_steps);
}


/* ---- cost, tangent directional derivative, adjoint gradient ---- */

static double cost_primal(double klf, double perc_limiter, int n_sub_setting)
{
    double cost;

    dsbm_cost_from_params(klf, perc_limiter, problem.theta_in, problem.n_steps,
                          n_sub_setting, problem.rain_mm_per_h, problem.pet_mm_per_h,
                          problem.observed_lateral_m, problem.observed_percolation_m,
                          problem.sst_lateral_m2, problem.sst_percolation_m2, &cost);
    return cost;
}

static double cost_tangent_direction(double klf, double perc_limiter, int n_sub_setting,
                                     double seed_klf, double seed_perc)
{
    double cost;
    double directional_derivative;

    cost_tangent_from_tapenade(klf, perc_limiter, seed_klf, seed_perc, problem.theta_in,
                               problem.n_steps, n_sub_setting,
                               problem.rain_mm_per_h, problem.pet_mm_per_h,
                               problem.observed_lateral_m, problem.observed_percolation_m,
                               problem.sst_lateral_m2, problem.sst_percolation_m2,
                               &cost, &directional_derivative);
    return directional_derivative;
}

static void cost_gradient_adjoint(double klf, double perc_limiter, int n_sub_setting,
                                  double *d_cost_d_klf, double *d_cost_d_perc)
{
    cost_adjoint_from_tapenade(klf, perc_limiter, problem.theta_in, problem.n_steps,
                               n_sub_setting, problem.rain_mm_per_h, problem.pet_mm_per_h,
                               problem.observed_lateral_m, problem.observed_percolation_m,
                               problem.sst_lateral_m2, problem.sst_percolation_m2,
                               1.0, d_cost_d_klf, d_cost_d_perc);
}

/* Gradient by the selected method.  Returns the cost (from the primal). */
static double cost_gradient_ad(double klf, double perc_limiter, int n_sub_setting,
                               double *d_cost_d_klf, double *d_cost_d_perc)
{
    if (gradient_method == 1) {
        cost_gradient_adjoint(klf, perc_limiter, n_sub_setting, d_cost_d_klf, d_cost_d_perc);
    } else {
        *d_cost_d_klf = cost_tangent_direction(klf, perc_limiter, n_sub_setting, 1.0, 0.0);
        *d_cost_d_perc = cost_tangent_direction(klf, perc_limiter, n_sub_setting, 0.0, 1.0);
    }
    return cost_primal(klf, perc_limiter, n_sub_setting);
}

static double relative_difference(double value, double reference)
{
    if (fabs(reference) > 1.0e-300) return (value - reference) / fabs(reference);
    return value - reference;
}



/* ---- E9a adjoint versus tangent, dot-product test ---- */

static void experiment_9a_adjoint_check(void)
{
    const double klf = 5.0e-4;
    const double perc_limiter = 0.8;
    const int n_sub_list[3] = {4, 12, 24};
    const double directions[4][2] = {{1.0, 0.0}, {0.0, 1.0}, {1.0e-3, 1.0}, {-2.0e-4, 0.7}};

    printf("====================================================================\n");
    printf("E9a  ADJOINT versus TANGENT at k_lf = %.1e m/h, perc_limiter = %.2f (1 year)\n",
           klf, perc_limiter);
    printf("====================================================================\n");

    for (int i_n = 0; i_n < 3; i_n++) {
        int n_sub = n_sub_list[i_n];
        double adjoint_klf;
        double adjoint_perc;
        double tangent_klf;
        double tangent_perc;

        cost_gradient_adjoint(klf, perc_limiter, n_sub, &adjoint_klf, &adjoint_perc);
        tangent_klf = cost_tangent_direction(klf, perc_limiter, n_sub, 1.0, 0.0);
        tangent_perc = cost_tangent_direction(klf, perc_limiter, n_sub, 0.0, 1.0);

        printf("\n--- fixed %d substeps ---\n", n_sub);
        printf("%-22s %24s %24s %10s\n", "", "adjoint", "tangent", "rel diff");
        printf("%-22s %+24.16e %+24.16e %10.2e\n", "dcost/dk_lf", adjoint_klf, tangent_klf,
               relative_difference(adjoint_klf, tangent_klf));
        printf("%-22s %+24.16e %+24.16e %10.2e\n", "dcost/dperc_limiter", adjoint_perc,
               tangent_perc, relative_difference(adjoint_perc, tangent_perc));

        printf("dot-product test, w = 1:  (J v)  versus  (J^T w) . v\n");
        for (int i_v = 0; i_v < 4; i_v++) {
            double v_klf = directions[i_v][0];
            double v_perc = directions[i_v][1];
            double tangent_value = cost_tangent_direction(klf, perc_limiter, n_sub, v_klf, v_perc);
            double adjoint_value = adjoint_klf * v_klf + adjoint_perc * v_perc;

            printf("   v = (%+.1e, %+.2f):  %+.16e  %+.16e  rel diff %.2e\n", v_klf, v_perc,
                   tangent_value, adjoint_value,
                   relative_difference(adjoint_value, tangent_value));
        }
    }
}


/* ---- E9b run time ---- */

static void experiment_9b_timing(void)
{
    const int n_repeats = 5;
    const double klf = 5.0e-4;
    const double perc_limiter = 0.8;
    const int n_sub = 4;
    double best_primal = 1.0e30;
    double best_tangent = 1.0e30;
    double best_adjoint = 1.0e30;
    double d1;
    double d2;

    printf("\n====================================================================\n");
    printf("E9b  RUN TIME, %d hours, fixed %d substeps (best of %d)\n",
           problem.n_steps, n_sub, n_repeats);
    printf("====================================================================\n");

    for (int i_r = 0; i_r < n_repeats; i_r++) {
        double t0;
        double elapsed;

        t0 = wall_clock_seconds();
        cost_primal(klf, perc_limiter, n_sub);
        elapsed = wall_clock_seconds() - t0;
        if (elapsed < best_primal) best_primal = elapsed;

        t0 = wall_clock_seconds();
        cost_tangent_direction(klf, perc_limiter, n_sub, 1.0, 0.0);
        elapsed = wall_clock_seconds() - t0;
        if (elapsed < best_tangent) best_tangent = elapsed;

        t0 = wall_clock_seconds();
        cost_gradient_adjoint(klf, perc_limiter, n_sub, &d1, &d2);
        elapsed = wall_clock_seconds() - t0;
        if (elapsed < best_adjoint) best_adjoint = elapsed;
    }

    printf("primal cost function         %9.4f s   (1.00 x primal)\n", best_primal);
    printf("tangent, one direction       %9.4f s   (%.2f x primal)\n", best_tangent,
           best_tangent / best_primal);
    printf("adjoint, whole gradient      %9.4f s   (%.2f x primal)\n", best_adjoint,
           best_adjoint / best_primal);
    printf("\nGradient of N parameters: tangent needs N runs, adjoint one.\n");
    printf("Break-even N (adjoint time / tangent time) = %.2f\n", best_adjoint / best_tangent);
    printf("   N = 2  : tangent %.4f s   adjoint %.4f s\n", 2.0 * best_tangent, best_adjoint);
    printf("   N = 12 : tangent %.4f s   adjoint %.4f s  (projected)\n", 12.0 * best_tangent,
           best_adjoint);
}


/* ---- E9c memory ---- */

static void measure_adjoint_memory(int n_hours)
{
    long rss_before;
    long rss_after;
    double d1;
    double d2;
    double t0;
    double elapsed;
    double substeps;

    set_up_problem(n_hours);
    rss_before = peak_resident_memory_kilobytes();
    t0 = wall_clock_seconds();
    cost_gradient_adjoint(5.0e-4, 0.8, 4, &d1, &d2);
    elapsed = wall_clock_seconds() - t0;
    rss_after = peak_resident_memory_kilobytes();
    substeps = 4.0 * (double)n_hours;

    printf("%8d hours: peak memory growth %10.1f MB  (%7.1f bytes per substep)   adjoint time %.3f s\n",
           n_hours, (double)(rss_after - rss_before) / 1024.0,
           1024.0 * (double)(rss_after - rss_before) / substeps, elapsed);
}

static void experiment_9c_memory(void)
{
    printf("\n====================================================================\n");
    printf("E9c  ADJOINT MEMORY (peak resident set growth during the first adjoint\n");
    printf("     run of each length; fixed 4 substeps)\n");
    printf("====================================================================\n");
    printf("peak resident memory before any adjoint run: %.1f MB\n",
           (double)rss_before_first_adjoint_kb / 1024.0);
    printf("%8d hours: peak memory growth %10.1f MB  (%7.1f bytes per substep)   adjoint time %.3f s\n",
           HOURS_ONE_YEAR,
           (double)(rss_after_first_adjoint_kb - rss_before_first_adjoint_kb) / 1024.0,
           1024.0 * (double)(rss_after_first_adjoint_kb - rss_before_first_adjoint_kb)
               / (4.0 * (double)HOURS_ONE_YEAR),
           first_adjoint_seconds);
    measure_adjoint_memory(n_hours_available);
    printf("(the full-record line is growth beyond the 1-year peak; its total tape is\n");
    printf(" the sum of both lines)\n");
    set_up_problem(HOURS_ONE_YEAR);
}


/* ---- E9d BFGS calibration ---- */

static void clamp_parameters(double *log_klf, double *perc_limiter)
{
    if (*log_klf < log(KLF_LOWER_BOUND_M_PER_H)) *log_klf = log(KLF_LOWER_BOUND_M_PER_H);
    if (*log_klf > log(KLF_UPPER_BOUND_M_PER_H)) *log_klf = log(KLF_UPPER_BOUND_M_PER_H);
    if (*perc_limiter < PERC_LIMITER_LOWER_BOUND) *perc_limiter = PERC_LIMITER_LOWER_BOUND;
    if (*perc_limiter > PERC_LIMITER_UPPER_BOUND) *perc_limiter = PERC_LIMITER_UPPER_BOUND;
}

/* cost and gradient in u = (ln k_lf, perc_limiter) */
static double cost_and_gradient_u(const double *u, int n_sub, double *gradient_u)
{
    double klf = exp(u[0]);
    double d_cost_d_klf;
    double d_cost_d_perc;
    double cost;

    cost = cost_gradient_ad(klf, u[1], n_sub, &d_cost_d_klf, &d_cost_d_perc);
    gradient_u[0] = klf * d_cost_d_klf;
    gradient_u[1] = d_cost_d_perc;
    return cost;
}

static void calibrate_bfgs(double klf_start, double perc_start, int n_sub)
{
    double start_seconds = wall_clock_seconds();
    const int max_iterations = 60;
    double u[2];
    double gradient[2];
    double inverse_hessian[2][2] = {{1.0, 0.0}, {0.0, 1.0}};
    double cost;
    int n_cost_evaluations = 0;

    u[0] = log(klf_start);
    u[1] = perc_start;
    cost = cost_and_gradient_u(u, n_sub, gradient);
    n_cost_evaluations++;

    printf("\n--- %s gradient: start k_lf = %.2e m/h, perc_limiter = %.2f, fixed %d substeps ---\n",
           gradient_method_name(), klf_start, perc_start, n_sub);
    printf("%4s %16s %14s %12s %12s %6s\n", "iter", "cost", "k_lf (m/h)",
           "perc_limiter", "|gradient|", "evals");
    printf("%4d %16.8e %14.6e %12.6f %12.3e %6d\n", 0, cost, exp(u[0]), u[1],
           sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]), n_cost_evaluations);

    for (int iteration = 1; iteration <= max_iterations; iteration++) {
        double direction[2];
        double u_new[2];
        double gradient_new[2];
        double cost_new = 0.0;
        double step_length = 1.0;
        double slope;
        int accepted = 0;

        direction[0] = -(inverse_hessian[0][0] * gradient[0] + inverse_hessian[0][1] * gradient[1]);
        direction[1] = -(inverse_hessian[1][0] * gradient[0] + inverse_hessian[1][1] * gradient[1]);
        slope = direction[0] * gradient[0] + direction[1] * gradient[1];
        if (slope >= 0.0) {
            /* not a descent direction: reset to steepest descent */
            inverse_hessian[0][0] = 1.0;
            inverse_hessian[0][1] = 0.0;
            inverse_hessian[1][0] = 0.0;
            inverse_hessian[1][1] = 1.0;
            direction[0] = -gradient[0];
            direction[1] = -gradient[1];
            slope = direction[0] * gradient[0] + direction[1] * gradient[1];
        }

        /* backtracking line search with Armijo condition */
        for (int i_ls = 0; i_ls < 40; i_ls++) {
            u_new[0] = u[0] + step_length * direction[0];
            u_new[1] = u[1] + step_length * direction[1];
            clamp_parameters(&u_new[0], &u_new[1]);
            cost_new = cost_primal(exp(u_new[0]), u_new[1], n_sub);
            n_cost_evaluations++;
            if (cost_new <= cost + 1.0e-4 * step_length * slope) {
                accepted = 1;
                break;
            }
            step_length = 0.5 * step_length;
        }
        if (!accepted) {
            printf("line search failed at iteration %d\n", iteration);
            break;
        }

        cost_new = cost_and_gradient_u(u_new, n_sub, gradient_new);
        n_cost_evaluations++;
        {
            /* stop when the accepted step no longer reduces the cost
               meaningfully (roundoff floor of the cost function) */
            double relative_decrease = (cost - cost_new) / cost;
            if (relative_decrease < 1.0e-12) {
                u[0] = u_new[0];
                u[1] = u_new[1];
                cost = cost_new;
                printf("%4d %16.8e %14.6e %12.6f %12.3e %6d  (cost decrease < 1e-12 relative: converged)\n",
                       iteration, cost, exp(u[0]), u[1],
                       sqrt(gradient_new[0] * gradient_new[0] + gradient_new[1] * gradient_new[1]),
                       n_cost_evaluations);
                break;
            }
        }

        /* BFGS update of the inverse Hessian */
        {
            double s[2];
            double y[2];
            double sy;

            s[0] = u_new[0] - u[0];
            s[1] = u_new[1] - u[1];
            y[0] = gradient_new[0] - gradient[0];
            y[1] = gradient_new[1] - gradient[1];
            sy = s[0] * y[0] + s[1] * y[1];
            if (sy > 1.0e-300) {
                double Hy[2];
                double yHy;
                double rho = 1.0 / sy;

                Hy[0] = inverse_hessian[0][0] * y[0] + inverse_hessian[0][1] * y[1];
                Hy[1] = inverse_hessian[1][0] * y[0] + inverse_hessian[1][1] * y[1];
                yHy = y[0] * Hy[0] + y[1] * Hy[1];
                for (int r = 0; r < 2; r++) {
                    for (int c = 0; c < 2; c++) {
                        inverse_hessian[r][c] = inverse_hessian[r][c]
                            - rho * (Hy[r] * s[c] + s[r] * Hy[c])
                            + (rho * rho * yHy + rho) * s[r] * s[c];
                    }
                }
            }
        }

        u[0] = u_new[0];
        u[1] = u_new[1];
        gradient[0] = gradient_new[0];
        gradient[1] = gradient_new[1];
        cost = cost_new;

        printf("%4d %16.8e %14.6e %12.6f %12.3e %6d\n", iteration, cost, exp(u[0]), u[1],
               sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]), n_cost_evaluations);

        if (cost < 1.0e-20) break;
        if (sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) < 1.0e-12) break;
    }
    printf("gradient method: %s   wall-clock time for the whole calibration: %.3f s\n",
           gradient_method_name(), wall_clock_seconds() - start_seconds);
    printf("result: k_lf = %.8e m/h (truth %.1e, rel err %.2e), perc_limiter = %.8f (truth %.2f, rel err %.2e)\n",
           exp(u[0]), KLF_TRUE_M_PER_H, relative_difference(exp(u[0]), KLF_TRUE_M_PER_H),
           u[1], PERC_LIMITER_TRUE, relative_difference(u[1], PERC_LIMITER_TRUE));
}

static void experiment_9d_calibration(void)
{
    printf("\n====================================================================\n");
    printf("E9d  BFGS CALIBRATION WITH ADJOINT versus TANGENT GRADIENTS (1 year)\n");
    printf("====================================================================\n");
    gradient_method = 1;
    calibrate_bfgs(5.0e-4, 0.80, 4);
    gradient_method = 0;
    calibrate_bfgs(5.0e-4, 0.80, 4);
}


/* ---- main ---- */

int main(int argc, char **argv)
{
    const char *output_dir = ".";
    const char *forcing_path = "forcing/rain_pet_example.csv";
    const char *experiments = "abcd";

    if (argc > 1) output_dir = argv[1];
    if (argc > 2) forcing_path = argv[2];
    if (argc > 3) experiments = argv[3];
    (void)output_dir;

    n_hours_available = read_forcing_csv(forcing_path, MAX_HOURS,
                                         all_rain_mm_per_h, all_pet_mm_per_h);
    if (n_hours_available < HOURS_ONE_YEAR) {
        fprintf(stderr, "ERROR: need at least %d forcing hours\n", HOURS_ONE_YEAR);
        return 1;
    }
    set_up_problem(HOURS_ONE_YEAR);

    {
        double d1;
        double d2;
        double cost_check = cost_primal(1.0e-3, 0.7, 4);
        double tangent_cost;
        double directional;

        cost_tangent_from_tapenade(1.0e-3, 0.7, 1.0, 0.0, problem.theta_in, problem.n_steps, 4,
                                   problem.rain_mm_per_h, problem.pet_mm_per_h,
                                   problem.observed_lateral_m, problem.observed_percolation_m,
                                   problem.sst_lateral_m2, problem.sst_percolation_m2,
                                   &tangent_cost, &directional);
        if (fabs(tangent_cost - cost_check) > 1.0e-12 * fabs(cost_check)) {
            printf("CONSISTENCY FAILURE: primal %.17e tangent %.17e\n", cost_check, tangent_cost);
            printf("libsoil.a is probably stale.  Rebuild with: make veryclean ; make\n");
            return 2;
        }
        rss_before_first_adjoint_kb = peak_resident_memory_kilobytes();
        first_adjoint_seconds = wall_clock_seconds();
        cost_gradient_adjoint(1.0e-3, 0.7, 4, &d1, &d2);
        first_adjoint_seconds = wall_clock_seconds() - first_adjoint_seconds;
        rss_after_first_adjoint_kb = peak_resident_memory_kilobytes();
        printf("Consistency check passed: primal and tangent cost values agree.\n");
        printf("Forcing hours available: %d; problem: first %d hours.\n",
               n_hours_available, problem.n_steps);
    }

    if (strchr(experiments, 'a') != NULL) experiment_9a_adjoint_check();
    if (strchr(experiments, 'b') != NULL) experiment_9b_timing();
    if (strchr(experiments, 'c') != NULL) experiment_9c_memory();
    if (strchr(experiments, 'd') != NULL) experiment_9d_calibration();
    return 0;
}
