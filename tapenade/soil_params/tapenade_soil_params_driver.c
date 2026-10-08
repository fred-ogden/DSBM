/*
 * tapenade_soil_params_driver.c
 *
 * E10: four calibration parameters, two of them soil hydraulic properties.
 *
 *   parameters p = (k_lf (m/h), perc_limiter (-), Ksat (cm/h), b (-))
 *
 * phi_sat is computed from Ksat inside the differentiated code with the
 * CFE3.1 regression phi_sat(cm) = 10.415 Ksat(cm/h)^(-0.3266), and theta_fc
 * (= theta_aet_eq_pet) from phi_sat and b, so both move with the
 * parameters.  See tapenade_soil_cost_primal.c.
 *
 * TWIN EXPERIMENT: the "observations" are DSBM's own hourly lateral flow
 * and percolation at the true parameters, so these experiments test the
 * derivative machinery (Tapenade tangent and adjoint through the
 * Ksat -> phi_sat -> theta_fc chain and the moving thresholds), not
 * whether the parameters are identifiable from discharge at a gauge.
 *
 *   E10a  Adjoint gradient against tangent gradient (four tangent runs)
 *         at the test point, fixed 4 and 12 substeps, plus the
 *         dot-product test w * (J v) = (J^T w) . v.
 *   E10b  Adjoint gradient against central finite differences for a
 *         range of step sizes (first time theta_fc moves with the
 *         parameters).
 *   E10c  Run time: primal, one tangent direction, full tangent gradient
 *         (four runs), adjoint (one run); adjoint memory for 1 year.
 *   E10d  BFGS calibration over u = (ln k_lf, perc_limiter, ln Ksat, b)
 *         with the adjoint gradient and with the tangent gradient.
 *         Ksat is treated as a logarithmic variable, as is usual in soil
 *         science, and k_lf as in E8/E9.
 *
 * Usage:
 *   tapenade_soil_params_test  output_directory  forcing_csv_file  [experiments]
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


/* ---- functions in tapenade_soil_cost_primal.c and the two adapters ---- */

double soil_cost_phi_sat_cm_from_ksat(double K_sat_cm_per_h);
double soil_cost_theta_fc(double theta_sat, double phi_sat_cm, double b_exp);

void dsbm_soil_cost_from_params(double klf_m_per_h,
                                double perc_limiter_0_to_1,
                                double K_sat_cm_per_h,
                                double b_exp,
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

void dsbm_soil_series_from_params(double klf_m_per_h,
                                  double perc_limiter_0_to_1,
                                  double K_sat_cm_per_h,
                                  double b_exp,
                                  const double *theta_in,
                                  int n_steps,
                                  int n_sub_setting,
                                  const double *rain_mm_per_h,
                                  const double *pet_mm_per_h,
                                  double *lateral_series_m,
                                  double *percolation_series_m);

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
                                     double *directional_derivative_of_cost);

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
                                     double gradient_of_cost[4]);


/* ---- experiment constants ---- */

#define N_PARAMETERS                4
#define MAX_HOURS                   80000
#define HOURS_ONE_YEAR              8760
#define N_SUB_FOR_OBSERVATIONS      4

#define SOIL_THETA_SAT              0.439
#define SOIL_DEPTH_M                2.0
#define INITIAL_STORAGE_M           0.7843498367003

/* true parameters: configs/soil_params.dat soil, E8 k_lf and perc_limiter */
#define KLF_TRUE_M_PER_H            2.0e-3
#define PERC_LIMITER_TRUE           0.5
#define K_SAT_TRUE_CM_PER_H         1.2168
#define B_EXP_TRUE                  4.05

/* calibration bounds */
#define KLF_LOWER_BOUND_M_PER_H     1.0e-5
#define KLF_UPPER_BOUND_M_PER_H     1.0
#define PERC_LIMITER_LOWER_BOUND    0.01
#define PERC_LIMITER_UPPER_BOUND    0.999
#define K_SAT_LOWER_BOUND_CM_PER_H  0.01
#define K_SAT_UPPER_BOUND_CM_PER_H  100.0
#define B_EXP_LOWER_BOUND           2.0
#define B_EXP_UPPER_BOUND           15.0

static const double disc_center_depth_m[NDISC] = {0.05, 0.25, 0.70, 1.50};
static const char *parameter_name[N_PARAMETERS] = {"k_lf", "perc_limiter", "Ksat", "b"};
static const double parameter_truth[N_PARAMETERS] = {
    KLF_TRUE_M_PER_H, PERC_LIMITER_TRUE, K_SAT_TRUE_CM_PER_H, B_EXP_TRUE};

/* test point for E10a-E10c, and the BFGS start of E10d */
static const double parameter_test_point[N_PARAMETERS] = {5.0e-4, 0.8, 3.0, 5.5};


/* ---- forcing (whole record) and the calibration problem ---- */

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

/* peak resident memory before and after the very first adjoint run */
static long rss_before_first_adjoint_kb = 0;
static long rss_after_first_adjoint_kb = 0;

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
    char *first_comma;
    double rain_value;
    double pet_value;
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
        first_comma = strchr(line, ',');
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
    double deviation;
    int i;

    for (i = 0; i < n; i++) mean_value = mean_value + series[i];
    mean_value = mean_value / (double)n;
    for (i = 0; i < n; i++) {
        deviation = series[i] - mean_value;
        sum_squares = sum_squares + deviation * deviation;
    }
    return sum_squares;
}

/*
 * Set up the synthetic-truth problem over the first n_hours of forcing.
 * The initial soil moisture is hydrostatic for the TRUE soil and is a
 * passive input: it does not change with the trial parameters.
 */
static void set_up_problem(int n_hours)
{
    double water_table_depth_m;

    problem.n_steps = n_hours;
    problem.rain_mm_per_h = all_rain_mm_per_h;
    problem.pet_mm_per_h = all_pet_mm_per_h;

    initialize_hydrostatic_from_storage(SOIL_DEPTH_M, SOIL_THETA_SAT,
                                        soil_cost_phi_sat_cm_from_ksat(K_SAT_TRUE_CM_PER_H),
                                        B_EXP_TRUE, disc_center_depth_m, INITIAL_STORAGE_M,
                                        &water_table_depth_m, problem.theta_in);

    dsbm_soil_series_from_params(KLF_TRUE_M_PER_H, PERC_LIMITER_TRUE,
                                 K_SAT_TRUE_CM_PER_H, B_EXP_TRUE, problem.theta_in,
                                 problem.n_steps, N_SUB_FOR_OBSERVATIONS,
                                 problem.rain_mm_per_h, problem.pet_mm_per_h,
                                 problem.observed_lateral_m, problem.observed_percolation_m);

    problem.sst_lateral_m2 = sum_of_squared_deviations(problem.observed_lateral_m, problem.n_steps);
    problem.sst_percolation_m2 = sum_of_squared_deviations(problem.observed_percolation_m, problem.n_steps);
}


/* ---- cost, tangent directional derivative, adjoint gradient ---- */

static double cost_primal(const double p[N_PARAMETERS], int n_sub_setting)
{
    double cost;

    dsbm_soil_cost_from_params(p[0], p[1], p[2], p[3], problem.theta_in, problem.n_steps,
                               n_sub_setting, problem.rain_mm_per_h, problem.pet_mm_per_h,
                               problem.observed_lateral_m, problem.observed_percolation_m,
                               problem.sst_lateral_m2, problem.sst_percolation_m2, &cost);
    return cost;
}

static double cost_tangent_direction(const double p[N_PARAMETERS], int n_sub_setting,
                                     const double seed[N_PARAMETERS])
{
    double cost;
    double directional_derivative;

    soil_cost_tangent_from_tapenade(p, seed, problem.theta_in, problem.n_steps, n_sub_setting,
                                    problem.rain_mm_per_h, problem.pet_mm_per_h,
                                    problem.observed_lateral_m, problem.observed_percolation_m,
                                    problem.sst_lateral_m2, problem.sst_percolation_m2,
                                    &cost, &directional_derivative);
    return directional_derivative;
}

static void cost_gradient_tangent(const double p[N_PARAMETERS], int n_sub_setting,
                                  double gradient[N_PARAMETERS])
{
    double seed[N_PARAMETERS];
    int i;
    int j;

    for (i = 0; i < N_PARAMETERS; i++) {
        for (j = 0; j < N_PARAMETERS; j++) seed[j] = 0.0;
        seed[i] = 1.0;
        gradient[i] = cost_tangent_direction(p, n_sub_setting, seed);
    }
}

static void cost_gradient_adjoint(const double p[N_PARAMETERS], int n_sub_setting,
                                  double gradient[N_PARAMETERS])
{
    soil_cost_adjoint_from_tapenade(p, problem.theta_in, problem.n_steps, n_sub_setting,
                                    problem.rain_mm_per_h, problem.pet_mm_per_h,
                                    problem.observed_lateral_m, problem.observed_percolation_m,
                                    problem.sst_lateral_m2, problem.sst_percolation_m2,
                                    1.0, gradient);
}

/* Gradient by the selected method.  Returns the cost (from the primal). */
static double cost_gradient_ad(const double p[N_PARAMETERS], int n_sub_setting,
                               double gradient[N_PARAMETERS])
{
    if (gradient_method == 1) {
        cost_gradient_adjoint(p, n_sub_setting, gradient);
    } else {
        cost_gradient_tangent(p, n_sub_setting, gradient);
    }
    return cost_primal(p, n_sub_setting);
}

static double relative_difference(double value, double reference)
{
    if (fabs(reference) > 1.0e-300) return (value - reference) / fabs(reference);
    return value - reference;
}

static void print_derived_soil(const char *label, const double p[N_PARAMETERS])
{
    double phi_sat_cm = soil_cost_phi_sat_cm_from_ksat(p[2]);
    double theta_fc = soil_cost_theta_fc(SOIL_THETA_SAT, phi_sat_cm, p[3]);

    printf("%s Ksat = %.4f cm/h, b = %.3f  ->  phi_sat = %.3f cm, theta_fc = %.4f\n",
           label, p[2], p[3], phi_sat_cm, theta_fc);
}


/* ---- E10a adjoint versus tangent, dot-product test ---- */

static void experiment_10a_adjoint_check(void)
{
    const int n_sub_list[2] = {4, 12};
    const double directions[4][N_PARAMETERS] = {
        {1.0, 0.0, 0.0, 0.0},
        {1.0e-3, 1.0, 0.1, 0.5},
        {-2.0e-4, 0.3, -0.5, 1.0},
        {0.0, 0.0, 1.0, -1.0}};
    double adjoint_gradient[N_PARAMETERS];
    double tangent_gradient[N_PARAMETERS];
    double tangent_value;
    double adjoint_value;
    int i_n;
    int i;
    int i_v;

    printf("====================================================================\n");
    printf("E10a ADJOINT versus TANGENT, 4 parameters (1 year)\n");
    printf("====================================================================\n");
    print_derived_soil("test point:", parameter_test_point);
    printf("            k_lf = %.1e m/h, perc_limiter = %.2f\n",
           parameter_test_point[0], parameter_test_point[1]);
    print_derived_soil("truth:     ", parameter_truth);

    for (i_n = 0; i_n < 2; i_n++) {
        cost_gradient_adjoint(parameter_test_point, n_sub_list[i_n], adjoint_gradient);
        cost_gradient_tangent(parameter_test_point, n_sub_list[i_n], tangent_gradient);

        printf("\n--- fixed %d substeps ---\n", n_sub_list[i_n]);
        printf("%-22s %24s %24s %10s\n", "", "adjoint", "tangent", "rel diff");
        for (i = 0; i < N_PARAMETERS; i++) {
            printf("dcost/d%-15s %+24.16e %+24.16e %10.2e\n", parameter_name[i],
                   adjoint_gradient[i], tangent_gradient[i],
                   relative_difference(adjoint_gradient[i], tangent_gradient[i]));
        }
        printf("dot-product test, w = 1:  (J v)  versus  (J^T w) . v\n");
        for (i_v = 0; i_v < 4; i_v++) {
            tangent_value = cost_tangent_direction(parameter_test_point, n_sub_list[i_n],
                                                   directions[i_v]);
            adjoint_value = 0.0;
            for (i = 0; i < N_PARAMETERS; i++) {
                adjoint_value = adjoint_value + adjoint_gradient[i] * directions[i_v][i];
            }
            printf("   v%d:  %+.16e  %+.16e  rel diff %.2e\n", i_v + 1, tangent_value,
                   adjoint_value, relative_difference(adjoint_value, tangent_value));
        }
    }
    printf("(v1 = k_lf only; v2, v3 mixed; v4 = Ksat up and b down together)\n");
}


/* ---- E10b adjoint versus central finite differences ---- */

static void experiment_10b_finite_differences(void)
{
    const double relative_steps[5] = {1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};
    double adjoint_gradient[N_PARAMETERS];
    double p_plus[N_PARAMETERS];
    double p_minus[N_PARAMETERS];
    double step;
    double finite_difference;
    int i;
    int j;
    int i_h;

    printf("\n====================================================================\n");
    printf("E10b ADJOINT versus CENTRAL FINITE DIFFERENCES, fixed 4 substeps (1 year)\n");
    printf("     step = relative step * parameter value\n");
    printf("====================================================================\n");

    cost_gradient_adjoint(parameter_test_point, 4, adjoint_gradient);

    for (i = 0; i < N_PARAMETERS; i++) {
        printf("\ndcost/d%s   adjoint %+.12e\n", parameter_name[i], adjoint_gradient[i]);
        printf("   %14s %22s %12s\n", "relative step", "finite difference", "rel diff");
        for (i_h = 0; i_h < 5; i_h++) {
            for (j = 0; j < N_PARAMETERS; j++) {
                p_plus[j] = parameter_test_point[j];
                p_minus[j] = parameter_test_point[j];
            }
            step = relative_steps[i_h] * parameter_test_point[i];
            p_plus[i] = parameter_test_point[i] + step;
            p_minus[i] = parameter_test_point[i] - step;
            finite_difference = (cost_primal(p_plus, 4) - cost_primal(p_minus, 4)) / (2.0 * step);
            printf("   %14.0e %+22.12e %12.2e\n", relative_steps[i_h], finite_difference,
                   relative_difference(finite_difference, adjoint_gradient[i]));
        }
    }
    printf("\nLarge steps: truncation error; small steps: roundoff and the kinks where\n");
    printf("a threshold (theta_fc, theta_wp, caps) switches inside the step.\n");
}


/* ---- E10c run time and memory ---- */

static void experiment_10c_timing(void)
{
    const int n_repeats = 5;
    const int n_sub = 4;
    double best_primal = 1.0e30;
    double best_tangent = 1.0e30;
    double best_adjoint = 1.0e30;
    double seed[N_PARAMETERS] = {1.0, 0.0, 0.0, 0.0};
    double gradient[N_PARAMETERS];
    double t0;
    double elapsed;
    int i_r;

    printf("\n====================================================================\n");
    printf("E10c RUN TIME, %d hours, fixed %d substeps (best of %d), and memory\n",
           problem.n_steps, n_sub, n_repeats);
    printf("====================================================================\n");

    for (i_r = 0; i_r < n_repeats; i_r++) {
        t0 = wall_clock_seconds();
        cost_primal(parameter_test_point, n_sub);
        elapsed = wall_clock_seconds() - t0;
        if (elapsed < best_primal) best_primal = elapsed;

        t0 = wall_clock_seconds();
        cost_tangent_direction(parameter_test_point, n_sub, seed);
        elapsed = wall_clock_seconds() - t0;
        if (elapsed < best_tangent) best_tangent = elapsed;

        t0 = wall_clock_seconds();
        cost_gradient_adjoint(parameter_test_point, n_sub, gradient);
        elapsed = wall_clock_seconds() - t0;
        if (elapsed < best_adjoint) best_adjoint = elapsed;
    }

    printf("primal cost function            %9.4f s   (1.00 x primal)\n", best_primal);
    printf("tangent, one direction          %9.4f s   (%.2f x primal)\n", best_tangent,
           best_tangent / best_primal);
    printf("tangent, whole gradient (4 runs)%9.4f s   (%.2f x primal)\n", 4.0 * best_tangent,
           4.0 * best_tangent / best_primal);
    printf("adjoint, whole gradient         %9.4f s   (%.2f x primal)\n", best_adjoint,
           best_adjoint / best_primal);
    printf("adjoint speed-up over tangent for 4 parameters: %.2f\n",
           4.0 * best_tangent / best_adjoint);
    printf("adjoint memory, first 1-year run: peak resident growth %.1f MB (%.1f bytes per substep)\n",
           (double)(rss_after_first_adjoint_kb - rss_before_first_adjoint_kb) / 1024.0,
           1024.0 * (double)(rss_after_first_adjoint_kb - rss_before_first_adjoint_kb)
               / (4.0 * (double)HOURS_ONE_YEAR));
}


/* ---- E10d BFGS calibration over u = (ln k_lf, perc_limiter, ln Ksat, b) ---- */

static void parameters_from_u(const double u[N_PARAMETERS], double p[N_PARAMETERS])
{
    p[0] = exp(u[0]);
    p[1] = u[1];
    p[2] = exp(u[2]);
    p[3] = u[3];
}

static void clamp_u(double u[N_PARAMETERS])
{
    if (u[0] < log(KLF_LOWER_BOUND_M_PER_H)) u[0] = log(KLF_LOWER_BOUND_M_PER_H);
    if (u[0] > log(KLF_UPPER_BOUND_M_PER_H)) u[0] = log(KLF_UPPER_BOUND_M_PER_H);
    if (u[1] < PERC_LIMITER_LOWER_BOUND) u[1] = PERC_LIMITER_LOWER_BOUND;
    if (u[1] > PERC_LIMITER_UPPER_BOUND) u[1] = PERC_LIMITER_UPPER_BOUND;
    if (u[2] < log(K_SAT_LOWER_BOUND_CM_PER_H)) u[2] = log(K_SAT_LOWER_BOUND_CM_PER_H);
    if (u[2] > log(K_SAT_UPPER_BOUND_CM_PER_H)) u[2] = log(K_SAT_UPPER_BOUND_CM_PER_H);
    if (u[3] < B_EXP_LOWER_BOUND) u[3] = B_EXP_LOWER_BOUND;
    if (u[3] > B_EXP_UPPER_BOUND) u[3] = B_EXP_UPPER_BOUND;
}

/* cost and gradient in u; chain rule d/d(ln x) = x d/dx for k_lf and Ksat */
static double cost_and_gradient_u(const double u[N_PARAMETERS], int n_sub,
                                  double gradient_u[N_PARAMETERS])
{
    double p[N_PARAMETERS];
    double gradient_p[N_PARAMETERS];
    double cost;

    parameters_from_u(u, p);
    cost = cost_gradient_ad(p, n_sub, gradient_p);
    gradient_u[0] = p[0] * gradient_p[0];
    gradient_u[1] = gradient_p[1];
    gradient_u[2] = p[2] * gradient_p[2];
    gradient_u[3] = gradient_p[3];
    return cost;
}

static double vector_norm(const double v[N_PARAMETERS])
{
    double sum = 0.0;
    int i;

    for (i = 0; i < N_PARAMETERS; i++) sum = sum + v[i] * v[i];
    return sqrt(sum);
}

static void print_iteration(int iteration, double cost, const double u[N_PARAMETERS],
                            const double gradient[N_PARAMETERS], int n_evaluations)
{
    double p[N_PARAMETERS];

    parameters_from_u(u, p);
    printf("%4d %15.8e %12.5e %9.5f %9.5f %8.4f %10.3e %5d\n", iteration, cost, p[0], p[1],
           p[2], p[3], vector_norm(gradient), n_evaluations);
}

static void calibrate_bfgs(const double p_start[N_PARAMETERS], int n_sub)
{
    const int max_iterations = 200;
    double start_seconds = wall_clock_seconds();
    double u[N_PARAMETERS];
    double u_new[N_PARAMETERS];
    double p[N_PARAMETERS];
    double gradient[N_PARAMETERS];
    double gradient_new[N_PARAMETERS];
    double direction[N_PARAMETERS];
    double inverse_hessian[N_PARAMETERS][N_PARAMETERS];
    double s[N_PARAMETERS];
    double y[N_PARAMETERS];
    double Hy[N_PARAMETERS];
    double cost;
    double cost_new = 0.0;
    double step_length;
    double slope;
    double sy;
    double yHy;
    double rho;
    double relative_decrease;
    int n_cost_evaluations = 0;
    int accepted;
    int iteration;
    int i_ls;
    int r;
    int c;

    for (r = 0; r < N_PARAMETERS; r++) {
        for (c = 0; c < N_PARAMETERS; c++) inverse_hessian[r][c] = 0.0;
        inverse_hessian[r][r] = 1.0;
    }
    u[0] = log(p_start[0]);
    u[1] = p_start[1];
    u[2] = log(p_start[2]);
    u[3] = p_start[3];
    cost = cost_and_gradient_u(u, n_sub, gradient);
    n_cost_evaluations++;

    printf("\n--- %s gradient, fixed %d substeps ---\n", gradient_method_name(), n_sub);
    printf("%4s %15s %12s %9s %9s %8s %10s %5s\n", "iter", "cost", "k_lf(m/h)", "perc_lim",
           "Ksat", "b", "|grad u|", "evals");
    print_iteration(0, cost, u, gradient, n_cost_evaluations);

    for (iteration = 1; iteration <= max_iterations; iteration++) {
        slope = 0.0;
        for (r = 0; r < N_PARAMETERS; r++) {
            direction[r] = 0.0;
            for (c = 0; c < N_PARAMETERS; c++) {
                direction[r] = direction[r] - inverse_hessian[r][c] * gradient[c];
            }
            slope = slope + direction[r] * gradient[r];
        }
        if (slope >= 0.0) {
            /* not a descent direction: reset to steepest descent */
            slope = 0.0;
            for (r = 0; r < N_PARAMETERS; r++) {
                for (c = 0; c < N_PARAMETERS; c++) inverse_hessian[r][c] = 0.0;
                inverse_hessian[r][r] = 1.0;
                direction[r] = -gradient[r];
                slope = slope + direction[r] * gradient[r];
            }
        }

        /* backtracking line search with Armijo condition */
        accepted = 0;
        step_length = 1.0;
        for (i_ls = 0; i_ls < 40; i_ls++) {
            for (r = 0; r < N_PARAMETERS; r++) u_new[r] = u[r] + step_length * direction[r];
            clamp_u(u_new);
            parameters_from_u(u_new, p);
            cost_new = cost_primal(p, n_sub);
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

        /* stop when the accepted step no longer reduces the cost
           meaningfully (roundoff floor of the cost function) */
        relative_decrease = (cost - cost_new) / cost;
        if (relative_decrease < 1.0e-12) {
            for (r = 0; r < N_PARAMETERS; r++) u[r] = u_new[r];
            cost = cost_new;
            print_iteration(iteration, cost, u, gradient_new, n_cost_evaluations);
            printf("     (cost decrease < 1e-12 relative: converged)\n");
            break;
        }

        /* BFGS update of the inverse Hessian */
        sy = 0.0;
        for (r = 0; r < N_PARAMETERS; r++) {
            s[r] = u_new[r] - u[r];
            y[r] = gradient_new[r] - gradient[r];
            sy = sy + s[r] * y[r];
        }
        if (sy > 1.0e-300) {
            rho = 1.0 / sy;
            yHy = 0.0;
            for (r = 0; r < N_PARAMETERS; r++) {
                Hy[r] = 0.0;
                for (c = 0; c < N_PARAMETERS; c++) Hy[r] = Hy[r] + inverse_hessian[r][c] * y[c];
                yHy = yHy + y[r] * Hy[r];
            }
            for (r = 0; r < N_PARAMETERS; r++) {
                for (c = 0; c < N_PARAMETERS; c++) {
                    inverse_hessian[r][c] = inverse_hessian[r][c]
                        - rho * (Hy[r] * s[c] + s[r] * Hy[c])
                        + (rho * rho * yHy + rho) * s[r] * s[c];
                }
            }
        }

        for (r = 0; r < N_PARAMETERS; r++) {
            u[r] = u_new[r];
            gradient[r] = gradient_new[r];
        }
        cost = cost_new;
        print_iteration(iteration, cost, u, gradient, n_cost_evaluations);

        if (cost < 1.0e-20) break;
        if (vector_norm(gradient) < 1.0e-12) break;
    }

    parameters_from_u(u, p);
    printf("gradient method: %s   wall-clock time for the whole calibration: %.3f s\n",
           gradient_method_name(), wall_clock_seconds() - start_seconds);
    for (r = 0; r < N_PARAMETERS; r++) {
        printf("   %-13s %.10e   truth %.6e   rel err %+.2e\n", parameter_name[r], p[r],
               parameter_truth[r], relative_difference(p[r], parameter_truth[r]));
    }
    print_derived_soil("   derived:", p);
}

static void experiment_10d_calibration(void)
{
    printf("\n====================================================================\n");
    printf("E10d BFGS CALIBRATION, 4 parameters, u = (ln k_lf, perc_limiter, ln Ksat, b)\n");
    printf("     twin experiment (1 year); start = E10a test point\n");
    printf("====================================================================\n");
    gradient_method = 1;
    calibrate_bfgs(parameter_test_point, 4);
    gradient_method = 0;
    calibrate_bfgs(parameter_test_point, 4);
}


/* ---- main ---- */

int main(int argc, char **argv)
{
    const char *output_dir = ".";
    const char *forcing_path = "forcing/rain_pet_example.csv";
    const char *experiments = "abcd";
    const double check_point[N_PARAMETERS] = {1.0e-3, 0.7, 2.0, 5.0};
    const double check_seed[N_PARAMETERS] = {1.0, 0.0, 0.0, 0.0};
    double gradient[N_PARAMETERS];
    double cost_check;
    double tangent_cost;
    double directional;

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

    /* the generated code and libsoil.a must come from the same sources */
    cost_check = cost_primal(check_point, 4);
    soil_cost_tangent_from_tapenade(check_point, check_seed, problem.theta_in, problem.n_steps, 4,
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
    cost_gradient_adjoint(check_point, 4, gradient);
    rss_after_first_adjoint_kb = peak_resident_memory_kilobytes();
    printf("Consistency check passed: primal and tangent cost values agree.\n");
    printf("Forcing hours available: %d; problem: first %d hours.\n",
           n_hours_available, problem.n_steps);
    printf("TWIN EXPERIMENT: tests the derivatives, not identifiability from discharge.\n");

    if (strchr(experiments, 'a') != NULL) experiment_10a_adjoint_check();
    if (strchr(experiments, 'b') != NULL) experiment_10b_finite_differences();
    if (strchr(experiments, 'c') != NULL) experiment_10c_timing();
    if (strchr(experiments, 'd') != NULL) experiment_10d_calibration();
    return 0;
}
