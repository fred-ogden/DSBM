/*
 * tapenade_klf_driver.c
 *
 * Test driver for the k_lf (lateral subsurface flow rate constant)
 * Tapenade tangent-mode experiments on the DSBM kernel.
 *
 * The differentiated function is dsbm_lateral_from_klf() in
 * tapenade_klf_primal.c.  Its Tapenade tangent is reached only through
 * klf_tangent_from_tapenade() in tapenade_klf_tangent_adapter.c, so this
 * file never needs to know the exact generated signature and can include
 * the production soil headers without colliding with Tapenade's
 * generated *_d.h headers.
 *
 * Experiments (run in this order; each writes a CSV to the output dir):
 *
 *   E1  One timestep, baseline entering state, no rain, no PET.
 *       AD versus centered finite differences (FD) for all four
 *       lateral_by_disc_m[], all four theta_out[], percolation, AET,
 *       and precipitation excess.  Analytic check where n_sub = 1.
 *       Volume-conservation check applied to the derivative itself.
 *
 *   E2  One timestep, storage-cap kink.  For n_sub = 1 the lateral
 *       removal from disc i is capped at (theta_i - theta_fc)*dz_i when
 *           k_lf >= k_cap_i = dz_i * (theta_sat - theta_fc) / dt_sub
 *       which is independent of theta.  Above k_cap_i the disc is
 *       drained exactly to theta_fc and d(lateral_i)/dk_lf = 0.
 *       AD versus forward, backward, and centered FD on both sides of
 *       each k_cap_i, plus a log-spaced sweep written to CSV.
 *
 *   E3  Multi-timestep drydown (no rain, constant PET) from the baseline
 *       state.  ET and redistribution carry discs across theta_fc, so
 *       the theta > theta_fc lateral-flow branch switches off at a time
 *       that depends on k_lf.  k_lf sweep; for each point AD versus
 *       FD, plus flags showing whether the +h and -h runs used a
 *       different adaptive n_sub sequence or a different theta_fc
 *       crossing hour.
 *
 *   E4  Multi-timestep run with observed forcing (first 2208 hours of
 *       forcing/rain_pet_example.csv).  First representative
 *       calibration derivative of the static k_lf parameter.
 *
 *   E5  Adaptive n_sub boundaries with the E4 observed forcing:
 *       a) dense k_lf sweep with production adaptive n_sub, comparing the
 *          actual change between neighbouring k_lf values with the change
 *          predicted by the AD derivative, separately for intervals with
 *          and without an n_sub switch;
 *       b) the worst switch located by bisection: the jump in every
 *          output, AD on each side, and centered FD straddling it;
 *       c) the same sweep with fixed n_sub = 12;
 *       d) adaptive versus fixed n_sub = 12 values and derivatives.
 *
 *   E6  Exact exponential lateral removal (ctrl.lateral_analytic = 1)
 *       versus forward Euler: one-timestep check against the analytic
 *       derivative (including k_lf above the Euler storage cap), the E5
 *       sweeps repeated with the exponential scheme, and how much values
 *       and derivatives depend on adaptive versus fixed n_sub = 12 under
 *       each scheme.
 *
 * Usage:
 *   tapenade_klf_derivative_test  output_directory  forcing_csv_file  [experiments]
 *
 *   experiments is a string of digits, default "123456"; e.g. "6" runs E6 only.
 *
 * Terminology: soil moisture; discs; cost function.
 * ASCII only.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "soil_data_types.h"
#include "soil_helpers.h"


/* ------------------------------------------------------------------ */
/* Functions defined in tapenade_klf_primal.c and the tangent adapter  */
/* ------------------------------------------------------------------ */

double klf_experiment_theta_fc(void);
double klf_experiment_phi_sat_cm(void);

void dsbm_lateral_from_klf(double klf_m_per_h,
                           const double *theta_in,
                           int n_steps,
                           int n_sub_fixed,
                           int lateral_analytic,
                           const double *rain_mm_per_h,
                           const double *pet_mm_per_h,
                           double *lateral_total_by_disc_m,
                           double *theta_out,
                           double *percolation_total_m,
                           double *aet_total_m,
                           double *precipitation_excess_total_m,
                           int *n_sub_used_by_step);

void klf_tangent_from_tapenade(double klf_m_per_h,
                               const double *theta_in,
                               int n_steps,
                               int n_sub_fixed,
                               int lateral_analytic,
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
                               int *n_sub_used_by_step);


/* ------------------------------------------------------------------ */
/* Experiment constants                                                */
/* ------------------------------------------------------------------ */

#define THETA_SAT_M3_PER_M3          0.439
#define SOIL_DEPTH_M                 2.0
#define BASELINE_INIT_STORAGE_M      0.7843498367003
#define BASELINE_KLF_M_PER_H         0.000005

/* output vector layout: 4 lateral, 4 theta_out, perc, AET, excess */
#define IDX_LATERAL_0     0
#define IDX_THETA_0       (NDISC)
#define IDX_PERC          (2 * NDISC)
#define IDX_AET           (2 * NDISC + 1)
#define IDX_EXCESS        (2 * NDISC + 2)
#define N_OUTPUTS         (2 * NDISC + 3)

#define E3_N_STEPS_H              480
#define E3_PET_MM_PER_H           0.30
#define E3_N_SWEEP                121
#define E4_N_STEPS_H              2208
#define MAX_STEPS                 (E4_N_STEPS_H)

static const double disc_thickness_m[NDISC] = {0.10, 0.30, 0.60, 1.00};
static const double disc_center_depth_m[NDISC] = {0.05, 0.25, 0.70, 1.50};

static const char *output_label[N_OUTPUTS] = {
    "lateral_by_disc_m[0]", "lateral_by_disc_m[1]",
    "lateral_by_disc_m[2]", "lateral_by_disc_m[3]",
    "theta_out[0]", "theta_out[1]", "theta_out[2]", "theta_out[3]",
    "percolation_to_gw_m", "AET_m", "precip_excess_m"
};

/* scale below which an AD value is treated as zero for relative errors */
#define TINY_DERIVATIVE_SCALE     1.0e-14


/* ------------------------------------------------------------------ */
/* Packing helpers                                                     */
/* ------------------------------------------------------------------ */

/*
 * Run the primal wrapper and pack all dependents into one vector.
 * n_sub_fixed = 0 uses production adaptive substeps.
 */
static void evaluate_primal_nsub(double klf_m_per_h,
                            const double *theta_in,
                            int n_steps,
                            int n_sub_fixed,
                            int lateral_analytic,
                            const double *rain_mm_per_h,
                            const double *pet_mm_per_h,
                            double *output_vector,
                            int *n_sub_used_by_step)
{
    double lateral_m[NDISC];
    double theta_out[NDISC];
    double perc_m;
    double aet_m;
    double excess_m;
    int i_disc;

    dsbm_lateral_from_klf(klf_m_per_h, theta_in, n_steps, n_sub_fixed, lateral_analytic,
                          rain_mm_per_h, pet_mm_per_h,
                          lateral_m, theta_out, &perc_m, &aet_m, &excess_m,
                          n_sub_used_by_step);

    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        output_vector[IDX_LATERAL_0 + i_disc] = lateral_m[i_disc];
        output_vector[IDX_THETA_0 + i_disc] = theta_out[i_disc];
    }
    output_vector[IDX_PERC] = perc_m;
    output_vector[IDX_AET] = aet_m;
    output_vector[IDX_EXCESS] = excess_m;
}

/* Run the Tapenade tangent (seed d(klf)=1) and pack values and derivatives. */
static void evaluate_tangent_nsub(double klf_m_per_h,
                             const double *theta_in,
                             int n_steps,
                             int n_sub_fixed,
                             int lateral_analytic,
                             const double *rain_mm_per_h,
                             const double *pet_mm_per_h,
                             double *output_vector,
                             double *d_output_d_klf,
                             int *n_sub_used_by_step)
{
    double lateral_m[NDISC];
    double d_lateral_m[NDISC];
    double theta_out[NDISC];
    double d_theta_out[NDISC];
    double perc_m;
    double d_perc_m;
    double aet_m;
    double d_aet_m;
    double excess_m;
    double d_excess_m;
    int i_disc;

    klf_tangent_from_tapenade(klf_m_per_h, theta_in, n_steps, n_sub_fixed, lateral_analytic,
                              rain_mm_per_h, pet_mm_per_h,
                              lateral_m, d_lateral_m,
                              theta_out, d_theta_out,
                              &perc_m, &d_perc_m,
                              &aet_m, &d_aet_m,
                              &excess_m, &d_excess_m,
                              n_sub_used_by_step);

    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        output_vector[IDX_LATERAL_0 + i_disc] = lateral_m[i_disc];
        output_vector[IDX_THETA_0 + i_disc] = theta_out[i_disc];
        d_output_d_klf[IDX_LATERAL_0 + i_disc] = d_lateral_m[i_disc];
        d_output_d_klf[IDX_THETA_0 + i_disc] = d_theta_out[i_disc];
    }
    output_vector[IDX_PERC] = perc_m;
    output_vector[IDX_AET] = aet_m;
    output_vector[IDX_EXCESS] = excess_m;
    d_output_d_klf[IDX_PERC] = d_perc_m;
    d_output_d_klf[IDX_AET] = d_aet_m;
    d_output_d_klf[IDX_EXCESS] = d_excess_m;
}

/* Production adaptive-substep versions used by E1-E4. */
static void evaluate_primal(double klf_m_per_h,
                            const double *theta_in,
                            int n_steps,
                            const double *rain_mm_per_h,
                            const double *pet_mm_per_h,
                            double *output_vector,
                            int *n_sub_used_by_step)
{
    evaluate_primal_nsub(klf_m_per_h, theta_in, n_steps, 0, 0,
                         rain_mm_per_h, pet_mm_per_h,
                         output_vector, n_sub_used_by_step);
}

static void evaluate_tangent(double klf_m_per_h,
                             const double *theta_in,
                             int n_steps,
                             const double *rain_mm_per_h,
                             const double *pet_mm_per_h,
                             double *output_vector,
                             double *d_output_d_klf,
                             int *n_sub_used_by_step)
{
    evaluate_tangent_nsub(klf_m_per_h, theta_in, n_steps, 0, 0,
                          rain_mm_per_h, pet_mm_per_h,
                          output_vector, d_output_d_klf, n_sub_used_by_step);
}

/*
 * Volume conservation applied to a derivative vector (m per (m/h)):
 *   sum dz*d(theta) + sum d(lateral) + d(perc) + d(AET) + d(excess)
 * Rain is fixed, so this must be zero if the tangent conserves volume.
 */
static double derivative_volume_residual(const double *d_output_d_klf)
{
    double residual;
    int i_disc;

    residual = 0.0;
    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        residual = residual
                 + disc_thickness_m[i_disc] * d_output_d_klf[IDX_THETA_0 + i_disc]
                 + d_output_d_klf[IDX_LATERAL_0 + i_disc];
    }
    residual = residual + d_output_d_klf[IDX_PERC]
                        + d_output_d_klf[IDX_AET]
                        + d_output_d_klf[IDX_EXCESS];
    return residual;
}

/* Relative difference of a FD estimate from AD, with absolute fallback. */
static double relative_difference_from_ad(double fd_value, double ad_value)
{
    double relative_difference;

    if (fabs(ad_value) > TINY_DERIVATIVE_SCALE) {
        relative_difference = (fd_value - ad_value) / fabs(ad_value);
    } else {
        relative_difference = fd_value - ad_value;
    }
    return relative_difference;
}

static const char *yes_no(int flag)
{
    if (flag) return "YES";
    return "no";
}

static const char *n_sub_note(int flag)
{
    if (flag) return "(n_sub changed)";
    return "";
}

/* 1 if two n_sub sequences differ anywhere, else 0 */
static int n_sub_sequences_differ(const int *a, const int *b, int n_steps)
{
    int i_step;

    for (i_step = 0; i_step < n_steps; i_step++) {
        if (a[i_step] != b[i_step]) return 1;
    }
    return 0;
}

/*
 * For each disc, the first hour (1-based) at whose end soil moisture is
 * at or below theta_fc; -1 if never.  Obtained by chaining one-step
 * primal calls, which is identical to the multi-step wrapper because the
 * wrapper recomputes psi and K from theta at the start of every step.
 */
static void first_fc_crossing_hours(double klf_m_per_h,
                                    const double *theta_in,
                                    int n_steps,
                                    const double *rain_mm_per_h,
                                    const double *pet_mm_per_h,
                                    int *first_crossing_hour_by_disc)
{
    double theta_now[NDISC];
    double output_vector[N_OUTPUTS];
    double theta_fc;
    int n_sub_one[1];
    int i_disc;
    int i_step;

    theta_fc = klf_experiment_theta_fc();
    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        theta_now[i_disc] = theta_in[i_disc];
        first_crossing_hour_by_disc[i_disc] = -1;
    }

    for (i_step = 0; i_step < n_steps; i_step++) {
        evaluate_primal(klf_m_per_h, theta_now, 1,
                        &rain_mm_per_h[i_step], &pet_mm_per_h[i_step],
                        output_vector, n_sub_one);
        for (i_disc = 0; i_disc < NDISC; i_disc++) {
            theta_now[i_disc] = output_vector[IDX_THETA_0 + i_disc];
            if (first_crossing_hour_by_disc[i_disc] < 0 &&
                theta_now[i_disc] <= theta_fc) {
                first_crossing_hour_by_disc[i_disc] = i_step + 1;
            }
        }
    }
}

static void baseline_entering_state(double *theta_in)
{
    double water_table_depth_m;

    initialize_hydrostatic_from_storage(SOIL_DEPTH_M,
                                        THETA_SAT_M3_PER_M3,
                                        klf_experiment_phi_sat_cm(),
                                        4.05,
                                        disc_center_depth_m,
                                        BASELINE_INIT_STORAGE_M,
                                        &water_table_depth_m,
                                        theta_in);
}

static FILE *open_csv(const char *output_dir, const char *file_name)
{
    char path[1024];
    FILE *fp;

    snprintf(path, sizeof(path), "%s/%s", output_dir, file_name);
    fp = fopen(path, "w");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", path);
        exit(1);
    }
    return fp;
}


/* ------------------------------------------------------------------ */
/* E1: one timestep, baseline state                                    */
/* ------------------------------------------------------------------ */

static void experiment_1_single_step(const char *output_dir)
{
    const double klf_values_m_per_h[2] = {BASELINE_KLF_M_PER_H, 1.0e-3};
    const double relative_eps[7] = {1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5,
                                    1.0e-6, 1.0e-7, 1.0e-8};
    double theta_in[NDISC];
    double rain_mm_per_h[1];
    double pet_mm_per_h[1];
    double y_ad[N_OUTPUTS];
    double dy_ad[N_OUTPUTS];
    double y_plus[N_OUTPUTS];
    double y_minus[N_OUTPUTS];
    double dy_fd[N_OUTPUTS];
    double theta_fc;
    int n_sub[1];
    int n_sub_pm[1];
    FILE *fp;

    theta_fc = klf_experiment_theta_fc();
    baseline_entering_state(theta_in);
    rain_mm_per_h[0] = 0.0;
    pet_mm_per_h[0] = 0.0;

    fp = open_csv(output_dir, "e1_single_step.csv");
    fprintf(fp, "klf_m_per_h,output,value,AD_derivative,rel_eps,FD_centered,rel_diff_FD_vs_AD\n");

    printf("====================================================================\n");
    printf("E1  ONE TIMESTEP, BASELINE ENTERING STATE, NO RAIN, NO PET\n");
    printf("====================================================================\n");
    printf("theta_sat = %.6f   theta_fc = %.6f\n", THETA_SAT_M3_PER_M3, theta_fc);
    for (int i_disc = 0; i_disc < NDISC; i_disc++) {
        printf("theta_in[%d] = %.6f  (theta_in - theta_fc = %+.6f)\n",
               i_disc, theta_in[i_disc], theta_in[i_disc] - theta_fc);
    }

    for (int i_klf = 0; i_klf < 2; i_klf++) {
        double klf = klf_values_m_per_h[i_klf];

        evaluate_tangent(klf, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                         y_ad, dy_ad, n_sub);

        printf("\n--- k_lf = %.6g m/h   n_sub = %d ---\n", klf, n_sub[0]);
        printf("%-22s %14s %16s %16s %16s %11s\n",
               "output", "value", "AD d/dk_lf", "analytic", "FD(1e-4 rel)", "FD-AD rel");

        /* FD at the reference perturbation */
        evaluate_primal(klf * (1.0 + 1.0e-4), theta_in, 1,
                        rain_mm_per_h, pet_mm_per_h, y_plus, n_sub_pm);
        evaluate_primal(klf * (1.0 - 1.0e-4), theta_in, 1,
                        rain_mm_per_h, pet_mm_per_h, y_minus, n_sub_pm);

        for (int k = 0; k < N_OUTPUTS; k++) {
            char analytic_text[32];

            dy_fd[k] = (y_plus[k] - y_minus[k]) / (2.0e-4 * klf);

            /*
             * Analytic value, only for n_sub = 1 and an uncapped disc:
             * lateral is exactly proportional to k_lf, so
             * d(lateral_i)/dk = lateral_i / k_lf,
             * d(theta_out_i)/dk = -(lateral_i / k_lf) / dz_i,
             * and percolation, AET, excess are independent of k_lf
             * because lateral removal is the last flux in the substep.
             */
            strcpy(analytic_text, "-");
            if (n_sub[0] == 1) {
                if (k < IDX_THETA_0) {
                    snprintf(analytic_text, sizeof(analytic_text), "%.9e",
                             y_ad[k] / klf);
                } else if (k < IDX_PERC) {
                    snprintf(analytic_text, sizeof(analytic_text), "%.9e",
                             -(y_ad[k - NDISC] / klf) /
                             disc_thickness_m[k - NDISC]);
                } else {
                    strcpy(analytic_text, "0");
                }
            }

            printf("%-22s %14.7e %16.9e %16s %16.9e %11.2e\n",
                   output_label[k], y_ad[k], dy_ad[k], analytic_text,
                   dy_fd[k], relative_difference_from_ad(dy_fd[k], dy_ad[k]));
        }

        printf("derivative volume residual (should be ~1e-16) = %.3e m/(m/h)\n",
               derivative_volume_residual(dy_ad));

        /* FD convergence: worst relative difference over all outputs */
        printf("FD convergence (worst |FD-AD|/|AD| over the 11 outputs):\n");
        for (int i_eps = 0; i_eps < 7; i_eps++) {
            double eps = relative_eps[i_eps];
            double worst = 0.0;
            const char *worst_label = "";

            evaluate_primal(klf * (1.0 + eps), theta_in, 1,
                            rain_mm_per_h, pet_mm_per_h, y_plus, n_sub_pm);
            evaluate_primal(klf * (1.0 - eps), theta_in, 1,
                            rain_mm_per_h, pet_mm_per_h, y_minus, n_sub_pm);
            for (int k = 0; k < N_OUTPUTS; k++) {
                double fd = (y_plus[k] - y_minus[k]) / (2.0 * eps * klf);
                double rd = relative_difference_from_ad(fd, dy_ad[k]);
                fprintf(fp, "%.9e,%s,%.17e,%.17e,%.1e,%.17e,%.6e\n",
                        klf, output_label[k], y_ad[k], dy_ad[k], eps, fd, rd);
                if (fabs(rd) > worst) {
                    worst = fabs(rd);
                    worst_label = output_label[k];
                }
            }
            printf("   rel_eps = %.0e   worst = %.3e  (%s)\n",
                   eps, worst, worst_label);
        }
    }
    fclose(fp);
}


/* ------------------------------------------------------------------ */
/* E2: one timestep, storage-cap kink                                  */
/* ------------------------------------------------------------------ */

static void experiment_2_storage_cap_kink(const char *output_dir)
{
    const double side_offsets[6] = {-1.0e-2, -1.0e-4, -1.0e-6,
                                    1.0e-6, 1.0e-4, 1.0e-2};
    const double fd_relative_step = 1.0e-8;
    double theta_in[NDISC];
    double rain_mm_per_h[1];
    double pet_mm_per_h[1];
    double y_ad[N_OUTPUTS];
    double dy_ad[N_OUTPUTS];
    double y_plus[N_OUTPUTS];
    double y_minus[N_OUTPUTS];
    double y_center[N_OUTPUTS];
    double theta_fc;
    double klf_cap_m_per_h[NDISC];
    int n_sub[1];
    int n_sub_pm[1];
    FILE *fp;

    theta_fc = klf_experiment_theta_fc();
    baseline_entering_state(theta_in);
    rain_mm_per_h[0] = 0.0;
    pet_mm_per_h[0] = 0.0;

    printf("\n====================================================================\n");
    printf("E2  ONE TIMESTEP, STORAGE-CAP KINK IN LATERAL REMOVAL\n");
    printf("====================================================================\n");
    printf("Predicted cap (n_sub = 1): k_cap_i = dz_i * (theta_sat - theta_fc) / dt_sub\n");
    for (int i_disc = 0; i_disc < NDISC; i_disc++) {
        klf_cap_m_per_h[i_disc] =
            disc_thickness_m[i_disc] * (THETA_SAT_M3_PER_M3 - theta_fc) / 1.0;
        printf("   disc %d   dz = %.2f m   k_cap = %.9e m/h\n",
               i_disc, disc_thickness_m[i_disc], klf_cap_m_per_h[i_disc]);
    }

    for (int i_disc = 0; i_disc < NDISC; i_disc++) {
        printf("\n--- disc %d near k_cap = %.6e m/h  (FD step = %.0e * k_lf) ---\n",
               i_disc, klf_cap_m_per_h[i_disc], fd_relative_step);
        printf("%12s %14s %15s %15s %15s %15s %6s\n",
               "k/k_cap - 1", "lateral_m", "AD", "FD forward", "FD backward",
               "FD centered", "n_sub");
        for (int i_off = 0; i_off < 6; i_off++) {
            double klf = klf_cap_m_per_h[i_disc] * (1.0 + side_offsets[i_off]);
            double h = klf * fd_relative_step;
            double fd_forward;
            double fd_backward;
            double fd_centered;

            evaluate_tangent(klf, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                             y_ad, dy_ad, n_sub);
            evaluate_primal(klf, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                            y_center, n_sub_pm);
            evaluate_primal(klf + h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                            y_plus, n_sub_pm);
            evaluate_primal(klf - h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                            y_minus, n_sub_pm);

            fd_forward = (y_plus[i_disc] - y_center[i_disc]) / h;
            fd_backward = (y_center[i_disc] - y_minus[i_disc]) / h;
            fd_centered = (y_plus[i_disc] - y_minus[i_disc]) / (2.0 * h);

            printf("%+12.0e %14.7e %15.8e %15.8e %15.8e %15.8e %6d\n",
                   side_offsets[i_off], y_ad[i_disc], dy_ad[i_disc],
                   fd_forward, fd_backward, fd_centered, n_sub[0]);
        }

        /* straddling centered FD exactly at the kink */
        {
            double klf = klf_cap_m_per_h[i_disc];
            double h = klf * 1.0e-3;
            double fd_straddle;

            evaluate_primal(klf + h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                            y_plus, n_sub_pm);
            evaluate_primal(klf - h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                            y_minus, n_sub_pm);
            fd_straddle = (y_plus[i_disc] - y_minus[i_disc]) / (2.0 * h);
            evaluate_tangent(klf, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                             y_ad, dy_ad, n_sub);
            printf("   exactly at k_cap: AD = %.8e   centered FD straddling (h=1e-3 k) = %.8e\n",
                   dy_ad[i_disc], fd_straddle);
        }
    }

    /* log-spaced sweep to CSV */
    fp = open_csv(output_dir, "e2_cap_sweep.csv");
    fprintf(fp, "klf_m_per_h");
    for (int i_disc = 0; i_disc < NDISC; i_disc++) {
        fprintf(fp, ",lateral_%d_m,AD_%d,FDc_%d", i_disc, i_disc, i_disc);
    }
    fprintf(fp, ",deriv_volume_residual\n");
    for (int i_pt = 0; i_pt <= 400; i_pt++) {
        double klf = pow(10.0, -4.0 + 4.0 * (double)i_pt / 400.0);
        double h = klf * 1.0e-6;

        evaluate_tangent(klf, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                         y_ad, dy_ad, n_sub);
        evaluate_primal(klf + h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                        y_plus, n_sub_pm);
        evaluate_primal(klf - h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                        y_minus, n_sub_pm);
        fprintf(fp, "%.9e", klf);
        for (int i_disc = 0; i_disc < NDISC; i_disc++) {
            fprintf(fp, ",%.12e,%.12e,%.12e", y_ad[i_disc], dy_ad[i_disc],
                    (y_plus[i_disc] - y_minus[i_disc]) / (2.0 * h));
        }
        fprintf(fp, ",%.3e\n", derivative_volume_residual(dy_ad));
    }
    fclose(fp);
    printf("\n(log-spaced sweep 1e-4 to 1 m/h written to e2_cap_sweep.csv)\n");
}


/* ------------------------------------------------------------------ */
/* Shared multi-step comparison of AD and FD at one k_lf value         */
/* ------------------------------------------------------------------ */

typedef struct {
    double value[N_OUTPUTS];
    double ad[N_OUTPUTS];
    double fd_centered[N_OUTPUTS];
    double fd_forward[N_OUTPUTS];
    double fd_backward[N_OUTPUTS];
    double ad_volume_residual;
    int n_sub_total;
    int n_sub_sequence_changed;     /* +h or -h run used different n_sub */
    int fc_crossing_hour_changed;   /* +h or -h run crossed theta_fc at a different hour */
} MultiStepComparison;

static void compare_multistep(double klf,
                              double fd_relative_step,
                              const double *theta_in,
                              int n_steps,
                              const double *rain_mm_per_h,
                              const double *pet_mm_per_h,
                              int check_fc_crossings,
                              MultiStepComparison *result)
{
    static int n_sub_center[MAX_STEPS];
    static int n_sub_plus[MAX_STEPS];
    static int n_sub_minus[MAX_STEPS];
    double y_plus[N_OUTPUTS];
    double y_minus[N_OUTPUTS];
    double h;
    int crossing_center[NDISC];
    int crossing_plus[NDISC];
    int crossing_minus[NDISC];

    h = klf * fd_relative_step;

    evaluate_tangent(klf, theta_in, n_steps, rain_mm_per_h, pet_mm_per_h,
                     result->value, result->ad, n_sub_center);
    evaluate_primal(klf + h, theta_in, n_steps, rain_mm_per_h, pet_mm_per_h,
                    y_plus, n_sub_plus);
    evaluate_primal(klf - h, theta_in, n_steps, rain_mm_per_h, pet_mm_per_h,
                    y_minus, n_sub_minus);

    for (int k = 0; k < N_OUTPUTS; k++) {
        result->fd_centered[k] = (y_plus[k] - y_minus[k]) / (2.0 * h);
        result->fd_forward[k] = (y_plus[k] - result->value[k]) / h;
        result->fd_backward[k] = (result->value[k] - y_minus[k]) / h;
    }
    result->ad_volume_residual = derivative_volume_residual(result->ad);

    result->n_sub_total = 0;
    for (int i_step = 0; i_step < n_steps; i_step++) {
        result->n_sub_total = result->n_sub_total + n_sub_center[i_step];
    }
    result->n_sub_sequence_changed = 0;
    if (n_sub_sequences_differ(n_sub_center, n_sub_plus, n_steps) ||
        n_sub_sequences_differ(n_sub_center, n_sub_minus, n_steps)) {
        result->n_sub_sequence_changed = 1;
    }

    result->fc_crossing_hour_changed = 0;
    if (check_fc_crossings) {
        first_fc_crossing_hours(klf, theta_in, n_steps, rain_mm_per_h,
                                pet_mm_per_h, crossing_center);
        first_fc_crossing_hours(klf + h, theta_in, n_steps, rain_mm_per_h,
                                pet_mm_per_h, crossing_plus);
        first_fc_crossing_hours(klf - h, theta_in, n_steps, rain_mm_per_h,
                                pet_mm_per_h, crossing_minus);
        for (int i_disc = 0; i_disc < NDISC; i_disc++) {
            if (crossing_plus[i_disc] != crossing_center[i_disc] ||
                crossing_minus[i_disc] != crossing_center[i_disc]) {
                result->fc_crossing_hour_changed = 1;
            }
        }
    }
}

static void print_multistep_table(double klf, const MultiStepComparison *r)
{
    printf("\n--- k_lf = %.6g m/h   total substeps = %d   n_sub changed under +-h: %s ---\n",
           klf, r->n_sub_total, yes_no(r->n_sub_sequence_changed));
    printf("%-22s %14s %16s %16s %11s\n",
           "output", "value", "AD d/dk_lf", "FD centered", "FD-AD rel");
    for (int k = 0; k < N_OUTPUTS; k++) {
        printf("%-22s %14.7e %16.9e %16.9e %11.2e\n",
               output_label[k], r->value[k], r->ad[k], r->fd_centered[k],
               relative_difference_from_ad(r->fd_centered[k], r->ad[k]));
    }
    printf("derivative volume residual = %.3e m/(m/h)\n", r->ad_volume_residual);
}


/* ------------------------------------------------------------------ */
/* E3: multi-timestep drydown, theta_fc threshold crossings            */
/* ------------------------------------------------------------------ */

static void experiment_3_drydown(const char *output_dir)
{
    static double rain_mm_per_h[E3_N_STEPS_H];
    static double pet_mm_per_h[E3_N_STEPS_H];
    double theta_in[NDISC];
    double theta_fc;
    int crossing_hours[NDISC];
    int n_points_tight = 0;
    int n_points_loose = 0;
    int n_points_bad = 0;
    double worst_rel = 0.0;
    double worst_klf = 0.0;
    MultiStepComparison r;
    FILE *fp;

    theta_fc = klf_experiment_theta_fc();
    baseline_entering_state(theta_in);
    for (int i_step = 0; i_step < E3_N_STEPS_H; i_step++) {
        rain_mm_per_h[i_step] = 0.0;
        pet_mm_per_h[i_step] = E3_PET_MM_PER_H;
    }

    printf("\n====================================================================\n");
    printf("E3  %d-HOUR DRYDOWN, NO RAIN, PET = %.2f mm/h, k_lf SWEEP\n",
           E3_N_STEPS_H, E3_PET_MM_PER_H);
    printf("====================================================================\n");
    printf("theta_fc = %.6f\n", theta_fc);

    {
        const double show_klf[3] = {BASELINE_KLF_M_PER_H, 1.0e-3, 1.0e-2};
        for (int i_show = 0; i_show < 3; i_show++) {
            first_fc_crossing_hours(show_klf[i_show], theta_in, E3_N_STEPS_H,
                                    rain_mm_per_h, pet_mm_per_h, crossing_hours);
            printf("k_lf = %.1e: first hour at/below theta_fc by disc = %d %d %d %d  (-1 = never)\n",
                   show_klf[i_show], crossing_hours[0], crossing_hours[1],
                   crossing_hours[2], crossing_hours[3]);
        }
        compare_multistep(1.0e-3, 1.0e-6, theta_in, E3_N_STEPS_H,
                          rain_mm_per_h, pet_mm_per_h, 1, &r);
        print_multistep_table(1.0e-3, &r);
    }

    fp = open_csv(output_dir, "e3_drydown_sweep.csv");
    fprintf(fp, "klf_m_per_h,lateral_total_m,AD_lateral_total,FDc_lateral_total,"
                "FDf_lateral_total,FDb_lateral_total,rel_diff,"
                "AD_theta0,FDc_theta0,AD_theta1,FDc_theta1,"
                "AD_perc,FDc_perc,AD_volume_residual,n_sub_total,"
                "n_sub_changed,fc_crossing_changed\n");

    printf("\nSweep of %d log-spaced k_lf values, 1e-4 to 1e-1 m/h, FD step 1e-6*k_lf.\n",
           E3_N_SWEEP);
    printf("Points where centered FD of total lateral flow disagrees with AD by > 1e-6:\n");
    printf("%14s %16s %16s %16s %16s %10s %5s %5s\n",
           "k_lf", "AD", "FD centered", "FD forward", "FD backward",
           "rel diff", "nsub", "fc_x");

    for (int i_pt = 0; i_pt < E3_N_SWEEP; i_pt++) {
        double klf = pow(10.0, -4.0 + 3.0 * (double)i_pt / (double)(E3_N_SWEEP - 1));
        double total_value = 0.0;
        double total_ad = 0.0;
        double total_fdc = 0.0;
        double total_fdf = 0.0;
        double total_fdb = 0.0;
        double rel;

        compare_multistep(klf, 1.0e-6, theta_in, E3_N_STEPS_H,
                          rain_mm_per_h, pet_mm_per_h, 1, &r);
        for (int i_disc = 0; i_disc < NDISC; i_disc++) {
            total_value = total_value + r.value[IDX_LATERAL_0 + i_disc];
            total_ad = total_ad + r.ad[IDX_LATERAL_0 + i_disc];
            total_fdc = total_fdc + r.fd_centered[IDX_LATERAL_0 + i_disc];
            total_fdf = total_fdf + r.fd_forward[IDX_LATERAL_0 + i_disc];
            total_fdb = total_fdb + r.fd_backward[IDX_LATERAL_0 + i_disc];
        }
        rel = relative_difference_from_ad(total_fdc, total_ad);

        if (fabs(rel) <= 1.0e-6) {
            n_points_tight++;
        } else if (fabs(rel) <= 1.0e-3) {
            n_points_loose++;
        } else {
            n_points_bad++;
        }
        if (fabs(rel) > worst_rel) {
            worst_rel = fabs(rel);
            worst_klf = klf;
        }
        if (fabs(rel) > 1.0e-6) {
            printf("%14.6e %16.9e %16.9e %16.9e %16.9e %10.2e %5s %5s\n",
                   klf, total_ad, total_fdc, total_fdf, total_fdb, rel,
                   yes_no(r.n_sub_sequence_changed),
                   yes_no(r.fc_crossing_hour_changed));
        }

        fprintf(fp, "%.9e,%.12e,%.12e,%.12e,%.12e,%.12e,%.6e,"
                    "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.3e,%d,%d,%d\n",
                klf, total_value, total_ad, total_fdc, total_fdf, total_fdb, rel,
                r.ad[IDX_THETA_0 + 0], r.fd_centered[IDX_THETA_0 + 0],
                r.ad[IDX_THETA_0 + 1], r.fd_centered[IDX_THETA_0 + 1],
                r.ad[IDX_PERC], r.fd_centered[IDX_PERC],
                r.ad_volume_residual, r.n_sub_total,
                r.n_sub_sequence_changed, r.fc_crossing_hour_changed);
    }
    fclose(fp);

    printf("Summary over %d points: |rel diff| <= 1e-6: %d   1e-6..1e-3: %d   > 1e-3: %d\n",
           E3_N_SWEEP, n_points_tight, n_points_loose, n_points_bad);
    printf("Worst |rel diff| = %.3e at k_lf = %.6e m/h\n", worst_rel, worst_klf);
    printf("(nsub = adaptive n_sub sequence differs between k-h, k, k+h;\n");
    printf(" fc_x = first theta_fc crossing hour of some disc differs between k-h, k, k+h)\n");
}


/* ------------------------------------------------------------------ */
/* E4: observed forcing, first representative calibration derivative   */
/* ------------------------------------------------------------------ */

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
    if (fgets(line, sizeof(line), fp) == NULL) {      /* header */
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

static void experiment_4_observed_forcing(const char *output_dir,
                                          const char *forcing_path)
{
    static double rain_mm_per_h[E4_N_STEPS_H];
    static double pet_mm_per_h[E4_N_STEPS_H];
    const double klf_values_m_per_h[3] = {BASELINE_KLF_M_PER_H, 1.0e-3, 1.0e-2};
    const double fd_steps[3] = {1.0e-4, 1.0e-6, 1.0e-8};
    double theta_in[NDISC];
    double rain_total_mm = 0.0;
    double pet_total_mm = 0.0;
    int n_steps;
    MultiStepComparison r;
    FILE *fp;

    n_steps = read_forcing_csv(forcing_path, E4_N_STEPS_H, rain_mm_per_h, pet_mm_per_h);
    if (n_steps <= 0) {
        printf("\nE4 skipped: could not read forcing file %s\n", forcing_path);
        return;
    }
    for (int i_step = 0; i_step < n_steps; i_step++) {
        rain_total_mm = rain_total_mm + rain_mm_per_h[i_step];
        pet_total_mm = pet_total_mm + pet_mm_per_h[i_step];
    }
    baseline_entering_state(theta_in);

    printf("\n====================================================================\n");
    printf("E4  OBSERVED FORCING, %d HOURS (rain %.1f mm, PET %.1f mm)\n",
           n_steps, rain_total_mm, pet_total_mm);
    printf("====================================================================\n");

    fp = open_csv(output_dir, "e4_observed_forcing.csv");
    fprintf(fp, "klf_m_per_h,fd_rel_step,output,value,AD,FD_centered,rel_diff,"
                "n_sub_total,n_sub_changed\n");

    for (int i_klf = 0; i_klf < 3; i_klf++) {
        double klf = klf_values_m_per_h[i_klf];

        compare_multistep(klf, 1.0e-6, theta_in, n_steps,
                          rain_mm_per_h, pet_mm_per_h, 0, &r);
        print_multistep_table(klf, &r);

        printf("worst |FD-AD|/|AD| over all outputs vs FD step:");
        for (int i_h = 0; i_h < 3; i_h++) {
            double worst = 0.0;

            compare_multistep(klf, fd_steps[i_h], theta_in, n_steps,
                              rain_mm_per_h, pet_mm_per_h, 0, &r);
            for (int k = 0; k < N_OUTPUTS; k++) {
                double rd = relative_difference_from_ad(r.fd_centered[k], r.ad[k]);
                if (fabs(rd) > worst) worst = fabs(rd);
                fprintf(fp, "%.9e,%.0e,%s,%.12e,%.12e,%.12e,%.6e,%d,%d\n",
                        klf, fd_steps[i_h], output_label[k], r.value[k], r.ad[k],
                        r.fd_centered[k], rd, r.n_sub_total,
                        r.n_sub_sequence_changed);
            }
            printf("  h=%.0e: %.2e%s", fd_steps[i_h], worst,
                   n_sub_note(r.n_sub_sequence_changed));
        }
        printf("\n");
    }
    fclose(fp);
}


/* ------------------------------------------------------------------ */
/* E5: adaptive n_sub boundaries versus fixed n_sub                    */
/* ------------------------------------------------------------------ */

#define E5_N_SWEEP              2000
#define E5_LOG10_KLF_MIN        (-4.0)
#define E5_LOG10_KLF_MAX        (-1.0)
#define E5_FIXED_N_SUB          12
#define E5_N_KEY_OUTPUTS        3

/*
 * Scalar outputs used to judge the sweep.  Total lateral flow is the sum
 * over discs; percolation and bottom-disc soil moisture are taken from
 * the output vector.
 */
static double total_lateral_from_vector(const double *output_vector)
{
    double total_m;
    int i_disc;

    total_m = 0.0;
    for (i_disc = 0; i_disc < NDISC; i_disc++) {
        total_m = total_m + output_vector[IDX_LATERAL_0 + i_disc];
    }
    return total_m;
}

static void key_outputs(const double *output_vector, double *key)
{
    key[0] = total_lateral_from_vector(output_vector);
    key[1] = output_vector[IDX_PERC];
    key[2] = output_vector[IDX_THETA_0 + NDISC - 1];
}

static const char *key_label[E5_N_KEY_OUTPUTS] = {
    "lateral_total_m", "percolation_to_gw_m", "theta_out[3]"
};

static int compare_doubles_ascending(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;

    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* number of timesteps whose n_sub differs between two sequences */
static int count_n_sub_differences(const int *a, const int *b, int n_steps)
{
    int n_differ;
    int i_step;

    n_differ = 0;
    for (i_step = 0; i_step < n_steps; i_step++) {
        if (a[i_step] != b[i_step]) n_differ++;
    }
    return n_differ;
}

/*
 * Summary of one k_lf sweep.  For each interval between neighbouring
 * k_lf values the actual change in a key output is compared with the
 * change predicted by trapezoidal integration of the AD derivative:
 *
 *   mismatch = [y(k2) - y(k1)] - 0.5*(dy/dk(k1) + dy/dk(k2))*(k2 - k1)
 *
 * and expressed as an EQUIVALENT RELATIVE CHANGE IN k_lf:
 *
 *   equiv_dk_over_k = |mismatch| / ( max|dy/dk| * k )
 *
 * i.e. the fractional change in k_lf that would move the output by the
 * same amount through the smooth derivative.  On smooth intervals this is
 * a tiny trapezoid error; at a slope kink it is at most of order the grid
 * spacing dk/k; at an n_sub switch it measures the jump.  A calibration
 * step smaller than this value cannot see the true gradient through the
 * jump.
 */
#define E5_TINY_DERIVATIVE      1.0e-12

typedef struct {
    int n_intervals;
    int n_intervals_with_n_sub_change;
    double max_rel_mismatch_smooth[E5_N_KEY_OUTPUTS];
    double max_rel_mismatch_switch[E5_N_KEY_OUTPUTS];
    double median_rel_mismatch_smooth[E5_N_KEY_OUTPUTS];
    double median_rel_mismatch_switch[E5_N_KEY_OUTPUTS];
    double worst_switch_klf_low;
    double worst_switch_klf_high;
    double worst_switch_rel_mismatch;
    long total_substeps_at_first_point;
} SweepSummary;

static void sweep_klf(int n_sub_fixed,
                      int lateral_analytic,
                      const double *theta_in,
                      int n_steps,
                      const double *rain_mm_per_h,
                      const double *pet_mm_per_h,
                      FILE *fp,
                      SweepSummary *summary)
{
    static int n_sub_previous[MAX_STEPS];
    static int n_sub_current[MAX_STEPS];
    static double smooth_values[E5_N_KEY_OUTPUTS][E5_N_SWEEP];
    static double switch_values[E5_N_KEY_OUTPUTS][E5_N_SWEEP];
    int n_smooth_values[E5_N_KEY_OUTPUTS];
    int n_switch_values[E5_N_KEY_OUTPUTS];
    double y_previous[N_OUTPUTS];
    double dy_previous[N_OUTPUTS];
    double y_current[N_OUTPUTS];
    double dy_current[N_OUTPUTS];
    double key_previous[E5_N_KEY_OUTPUTS];
    double key_current[E5_N_KEY_OUTPUTS];
    double dkey_previous[E5_N_KEY_OUTPUTS];
    double dkey_current[E5_N_KEY_OUTPUTS];
    double klf_previous;
    int i_pt;
    int k;

    summary->n_intervals = 0;
    summary->n_intervals_with_n_sub_change = 0;
    for (k = 0; k < E5_N_KEY_OUTPUTS; k++) {
        summary->max_rel_mismatch_smooth[k] = 0.0;
        summary->max_rel_mismatch_switch[k] = 0.0;
        summary->median_rel_mismatch_smooth[k] = 0.0;
        summary->median_rel_mismatch_switch[k] = 0.0;
        n_smooth_values[k] = 0;
        n_switch_values[k] = 0;
    }
    summary->worst_switch_klf_low = 0.0;
    summary->worst_switch_klf_high = 0.0;
    summary->worst_switch_rel_mismatch = -1.0;
    summary->total_substeps_at_first_point = 0;

    fprintf(fp, "klf_m_per_h,lateral_total_m,AD_lateral_total,percolation_m,AD_percolation,"
                "theta3,AD_theta3,n_sub_total,n_steps_changed_from_previous,"
                "equiv_dk_over_k_lateral,equiv_dk_over_k_percolation,equiv_dk_over_k_theta3\n");

    klf_previous = 0.0;
    for (i_pt = 0; i_pt < E5_N_SWEEP; i_pt++) {
        double log10_klf = E5_LOG10_KLF_MIN +
            (E5_LOG10_KLF_MAX - E5_LOG10_KLF_MIN) * (double)i_pt / (double)(E5_N_SWEEP - 1);
        double klf = pow(10.0, log10_klf);
        double rel_mismatch[E5_N_KEY_OUTPUTS];
        long n_sub_total;
        int n_steps_changed;
        int i_step;

        evaluate_tangent_nsub(klf, theta_in, n_steps, n_sub_fixed, lateral_analytic,
                              rain_mm_per_h, pet_mm_per_h,
                              y_current, dy_current, n_sub_current);
        key_outputs(y_current, key_current);
        dkey_current[0] = total_lateral_from_vector(dy_current);
        dkey_current[1] = dy_current[IDX_PERC];
        dkey_current[2] = dy_current[IDX_THETA_0 + NDISC - 1];

        n_sub_total = 0;
        for (i_step = 0; i_step < n_steps; i_step++) {
            n_sub_total = n_sub_total + n_sub_current[i_step];
        }
        if (i_pt == 0) summary->total_substeps_at_first_point = n_sub_total;

        n_steps_changed = 0;
        for (k = 0; k < E5_N_KEY_OUTPUTS; k++) rel_mismatch[k] = 0.0;

        if (i_pt > 0) {
            double dk = klf - klf_previous;
            double interval_worst = 0.0;

            n_steps_changed = count_n_sub_differences(n_sub_previous, n_sub_current, n_steps);
            summary->n_intervals++;
            if (n_steps_changed > 0) summary->n_intervals_with_n_sub_change++;

            for (k = 0; k < E5_N_KEY_OUTPUTS; k++) {
                double predicted = 0.5 * (dkey_previous[k] + dkey_current[k]) * dk;
                double actual = key_current[k] - key_previous[k];
                double derivative_scale = fabs(dkey_previous[k]);

                if (fabs(dkey_current[k]) > derivative_scale)
                    derivative_scale = fabs(dkey_current[k]);
                if (derivative_scale < E5_TINY_DERIVATIVE) continue;

                rel_mismatch[k] = fabs(actual - predicted) / (derivative_scale * klf);

                if (n_steps_changed > 0) {
                    if (rel_mismatch[k] > summary->max_rel_mismatch_switch[k])
                        summary->max_rel_mismatch_switch[k] = rel_mismatch[k];
                    switch_values[k][n_switch_values[k]] = rel_mismatch[k];
                    n_switch_values[k]++;
                } else {
                    if (rel_mismatch[k] > summary->max_rel_mismatch_smooth[k])
                        summary->max_rel_mismatch_smooth[k] = rel_mismatch[k];
                    smooth_values[k][n_smooth_values[k]] = rel_mismatch[k];
                    n_smooth_values[k]++;
                }
                if (rel_mismatch[k] > interval_worst) interval_worst = rel_mismatch[k];
            }

            if (n_steps_changed > 0 &&
                interval_worst > summary->worst_switch_rel_mismatch) {
                summary->worst_switch_rel_mismatch = interval_worst;
                summary->worst_switch_klf_low = klf_previous;
                summary->worst_switch_klf_high = klf;
            }
        }

        fprintf(fp, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%ld,%d,%.4e,%.4e,%.4e\n",
                klf, key_current[0], dkey_current[0], key_current[1], dkey_current[1],
                key_current[2], dkey_current[2], n_sub_total, n_steps_changed,
                rel_mismatch[0], rel_mismatch[1], rel_mismatch[2]);

        /* current becomes previous */
        klf_previous = klf;
        for (k = 0; k < N_OUTPUTS; k++) {
            y_previous[k] = y_current[k];
            dy_previous[k] = dy_current[k];
        }
        for (k = 0; k < E5_N_KEY_OUTPUTS; k++) {
            key_previous[k] = key_current[k];
            dkey_previous[k] = dkey_current[k];
        }
        for (i_step = 0; i_step < n_steps; i_step++) {
            n_sub_previous[i_step] = n_sub_current[i_step];
        }
    }
    (void)y_previous;
    (void)dy_previous;

    for (k = 0; k < E5_N_KEY_OUTPUTS; k++) {
        if (n_smooth_values[k] > 0) {
            qsort(smooth_values[k], (size_t)n_smooth_values[k], sizeof(double),
                  compare_doubles_ascending);
            summary->median_rel_mismatch_smooth[k] = smooth_values[k][n_smooth_values[k] / 2];
        }
        if (n_switch_values[k] > 0) {
            qsort(switch_values[k], (size_t)n_switch_values[k], sizeof(double),
                  compare_doubles_ascending);
            summary->median_rel_mismatch_switch[k] = switch_values[k][n_switch_values[k] / 2];
        }
    }
}

static void print_sweep_summary(const char *label, const SweepSummary *s)
{
    printf("\n%s: %d intervals, %d contain an n_sub switch in at least one timestep\n",
           label, s->n_intervals, s->n_intervals_with_n_sub_change);
    printf("   substeps per run at k_lf = 1e-4: %ld\n", s->total_substeps_at_first_point);
    printf("   equivalent dk/k of (actual - AD-predicted change) per interval:\n");
    printf("   %-22s %13s %13s %13s %13s\n", "output",
           "smooth median", "smooth max", "switch median", "switch max");
    for (int k = 0; k < E5_N_KEY_OUTPUTS; k++) {
        printf("   %-22s %13.3e %13.3e %13.3e %13.3e\n", key_label[k],
               s->median_rel_mismatch_smooth[k], s->max_rel_mismatch_smooth[k],
               s->median_rel_mismatch_switch[k], s->max_rel_mismatch_switch[k]);
    }
}

/*
 * Bisect inside [klf_low, klf_high] for the k_lf at which the n_sub
 * sequence switches, to near machine resolution.  Returns the two
 * adjacent k_lf values bracketing the switch.
 */
static void bisect_n_sub_switch(double klf_low, double klf_high,
                                int lateral_analytic,
                                const double *theta_in, int n_steps,
                                const double *rain_mm_per_h,
                                const double *pet_mm_per_h,
                                double *klf_minus, double *klf_plus)
{
    static int n_sub_low[MAX_STEPS];
    static int n_sub_mid[MAX_STEPS];
    double y_scratch[N_OUTPUTS];
    int iteration;

    evaluate_primal_nsub(klf_low, theta_in, n_steps, 0, lateral_analytic, rain_mm_per_h,
                         pet_mm_per_h, y_scratch, n_sub_low);

    for (iteration = 0; iteration < 200; iteration++) {
        double klf_mid = 0.5 * (klf_low + klf_high);

        if (klf_mid <= klf_low || klf_mid >= klf_high) break;
        evaluate_primal_nsub(klf_mid, theta_in, n_steps, 0, lateral_analytic, rain_mm_per_h,
                             pet_mm_per_h, y_scratch, n_sub_mid);
        if (count_n_sub_differences(n_sub_low, n_sub_mid, n_steps) == 0) {
            klf_low = klf_mid;
        } else {
            klf_high = klf_mid;
        }
    }
    *klf_minus = klf_low;
    *klf_plus = klf_high;
}

static void experiment_5_n_sub_boundaries(const char *output_dir,
                                          const char *forcing_path)
{
    static double rain_mm_per_h[E4_N_STEPS_H];
    static double pet_mm_per_h[E4_N_STEPS_H];
    static int n_sub_minus[MAX_STEPS];
    static int n_sub_plus[MAX_STEPS];
    const double compare_klf[3] = {BASELINE_KLF_M_PER_H, 1.0e-3, 1.0e-2};
    double theta_in[NDISC];
    double y_minus[N_OUTPUTS];
    double dy_minus[N_OUTPUTS];
    double y_plus[N_OUTPUTS];
    double dy_plus[N_OUTPUTS];
    double klf_minus;
    double klf_plus;
    int n_steps;
    SweepSummary adaptive_summary;
    SweepSummary fixed_summary;
    FILE *fp;

    n_steps = read_forcing_csv(forcing_path, E4_N_STEPS_H, rain_mm_per_h, pet_mm_per_h);
    if (n_steps <= 0) {
        printf("\nE5 skipped: could not read forcing file %s\n", forcing_path);
        return;
    }
    baseline_entering_state(theta_in);

    printf("\n====================================================================\n");
    printf("E5  ADAPTIVE n_sub BOUNDARIES, OBSERVED FORCING, %d HOURS\n", n_steps);
    printf("====================================================================\n");
    printf("Sweep of %d log-spaced k_lf values, 1e-4 to 1e-1 m/h (dk/k = %.2e).\n",
           E5_N_SWEEP,
           pow(10.0, (E5_LOG10_KLF_MAX - E5_LOG10_KLF_MIN) / (double)(E5_N_SWEEP - 1)) - 1.0);

    /* --- E5a adaptive sweep --- */
    fp = open_csv(output_dir, "e5_sweep_adaptive.csv");
    sweep_klf(0, 0, theta_in, n_steps, rain_mm_per_h, pet_mm_per_h, fp, &adaptive_summary);
    fclose(fp);
    print_sweep_summary("E5a ADAPTIVE n_sub", &adaptive_summary);

    /* --- E5b anatomy of the worst switch --- */
    if (adaptive_summary.n_intervals_with_n_sub_change > 0) {
        printf("\nE5b  WORST n_sub SWITCH, located by bisection\n");
        bisect_n_sub_switch(adaptive_summary.worst_switch_klf_low,
                            adaptive_summary.worst_switch_klf_high,
                            0,
                            theta_in, n_steps, rain_mm_per_h, pet_mm_per_h,
                            &klf_minus, &klf_plus);
        evaluate_tangent_nsub(klf_minus, theta_in, n_steps, 0, 0, rain_mm_per_h,
                              pet_mm_per_h, y_minus, dy_minus, n_sub_minus);
        evaluate_tangent_nsub(klf_plus, theta_in, n_steps, 0, 0, rain_mm_per_h,
                              pet_mm_per_h, y_plus, dy_plus, n_sub_plus);

        printf("switch between k_lf = %.17e\n", klf_minus);
        printf("           and k_lf = %.17e   (relative gap %.1e)\n",
               klf_plus, (klf_plus - klf_minus) / klf_minus);
        for (int i_step = 0; i_step < n_steps; i_step++) {
            if (n_sub_minus[i_step] != n_sub_plus[i_step]) {
                printf("   hour %d: n_sub %d -> %d   (rain %.3f mm/h)\n",
                       i_step + 1, n_sub_minus[i_step], n_sub_plus[i_step],
                       rain_mm_per_h[i_step]);
            }
        }
        printf("\n%-22s %15s %15s %15s %15s\n", "output", "jump J",
               "AD left", "AD right", "J/(AD*k)");
        for (int k = 0; k < N_OUTPUTS; k++) {
            double jump = y_plus[k] - y_minus[k];
            double equivalent_relative_dk = 0.0;
            char equivalent_text[32];

            strcpy(equivalent_text, "-");
            if (fabs(dy_minus[k]) > TINY_DERIVATIVE_SCALE) {
                equivalent_relative_dk = jump / (dy_minus[k] * klf_minus);
                snprintf(equivalent_text, sizeof(equivalent_text), "%.3e",
                         equivalent_relative_dk);
            }
            printf("%-22s %15.6e %15.8e %15.8e %15s\n", output_label[k],
                   jump, dy_minus[k], dy_plus[k], equivalent_text);
        }
        printf("(J/(AD*k) = the relative change in k_lf that would produce the same\n");
        printf(" change in that output through the smooth AD derivative)\n");

        printf("\nCentered FD of total lateral flow straddling the switch:\n");
        printf("%10s %18s %18s\n", "h/k", "FD centered", "AD (left side)");
        {
            const double straddle_h[6] = {1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12};
            double klf_switch = 0.5 * (klf_minus + klf_plus);
            static int n_sub_scratch[MAX_STEPS];

            for (int i_h = 0; i_h < 6; i_h++) {
                double h = klf_switch * straddle_h[i_h];
                double yp[N_OUTPUTS];
                double ym[N_OUTPUTS];

                evaluate_primal_nsub(klf_switch + h, theta_in, n_steps, 0, 0,
                                     rain_mm_per_h, pet_mm_per_h, yp, n_sub_scratch);
                evaluate_primal_nsub(klf_switch - h, theta_in, n_steps, 0, 0,
                                     rain_mm_per_h, pet_mm_per_h, ym, n_sub_scratch);
                printf("%10.0e %18.9e %18.9e\n", straddle_h[i_h],
                       (total_lateral_from_vector(yp) - total_lateral_from_vector(ym)) / (2.0 * h),
                       total_lateral_from_vector(dy_minus));
            }
        }
    }

    /* --- E5c fixed n_sub sweep --- */
    fp = open_csv(output_dir, "e5_sweep_fixed12.csv");
    sweep_klf(E5_FIXED_N_SUB, 0, theta_in, n_steps, rain_mm_per_h, pet_mm_per_h,
              fp, &fixed_summary);
    fclose(fp);
    print_sweep_summary("E5c FIXED n_sub = 12", &fixed_summary);

    /* --- E5d adaptive versus fixed at the E4 k_lf values --- */
    printf("\nE5d  ADAPTIVE versus FIXED n_sub = 12 at the E4 k_lf values\n");
    printf("%-10s %-22s %15s %15s %10s %15s %15s %10s\n", "k_lf", "output",
           "adaptive", "fixed12", "rel diff", "AD adaptive", "AD fixed12", "rel diff");
    for (int i_klf = 0; i_klf < 3; i_klf++) {
        double klf = compare_klf[i_klf];
        const int show_index[5] = {IDX_LATERAL_0, IDX_LATERAL_0 + NDISC - 1,
                                   IDX_THETA_0, IDX_THETA_0 + NDISC - 1, IDX_PERC};

        evaluate_tangent_nsub(klf, theta_in, n_steps, 0, 0, rain_mm_per_h,
                              pet_mm_per_h, y_minus, dy_minus, n_sub_minus);
        evaluate_tangent_nsub(klf, theta_in, n_steps, E5_FIXED_N_SUB, 0, rain_mm_per_h,
                              pet_mm_per_h, y_plus, dy_plus, n_sub_plus);
        for (int j = 0; j < 5; j++) {
            int k = show_index[j];
            printf("%-10.1e %-22s %15.8e %15.8e %10.2e %15.8e %15.8e %10.2e\n",
                   klf, output_label[k], y_minus[k], y_plus[k],
                   relative_difference_from_ad(y_plus[k], y_minus[k]),
                   dy_minus[k], dy_plus[k],
                   relative_difference_from_ad(dy_plus[k], dy_minus[k]));
        }
        {
            double lat_adaptive = total_lateral_from_vector(y_minus);
            double lat_fixed = total_lateral_from_vector(y_plus);
            double dlat_adaptive = total_lateral_from_vector(dy_minus);
            double dlat_fixed = total_lateral_from_vector(dy_plus);
            printf("%-10.1e %-22s %15.8e %15.8e %10.2e %15.8e %15.8e %10.2e\n",
                   klf, "lateral_total_m", lat_adaptive, lat_fixed,
                   relative_difference_from_ad(lat_fixed, lat_adaptive),
                   dlat_adaptive, dlat_fixed,
                   relative_difference_from_ad(dlat_fixed, dlat_adaptive));
            printf("%-10.1e derivative volume residual: adaptive %.2e   fixed12 %.2e\n",
                   klf, derivative_volume_residual(dy_minus),
                   derivative_volume_residual(dy_plus));
        }
    }
}


/* ------------------------------------------------------------------ */
/* E6: exact exponential lateral removal versus forward Euler          */
/* ------------------------------------------------------------------ */

/* jump across a located n_sub switch, as equivalent dk/k, for key outputs */
static void print_switch_jumps(const char *label, int lateral_analytic,
                               const SweepSummary *summary,
                               const double *theta_in, int n_steps,
                               const double *rain_mm_per_h,
                               const double *pet_mm_per_h)
{
    static int n_sub_minus[MAX_STEPS];
    static int n_sub_plus[MAX_STEPS];
    double y_minus[N_OUTPUTS];
    double dy_minus[N_OUTPUTS];
    double y_plus[N_OUTPUTS];
    double dy_plus[N_OUTPUTS];
    double klf_minus;
    double klf_plus;
    const int show_index[5] = {IDX_LATERAL_0, IDX_THETA_0, IDX_THETA_0 + 1,
                               IDX_THETA_0 + NDISC - 1, IDX_PERC};

    if (summary->n_intervals_with_n_sub_change == 0) {
        printf("%s: no n_sub switch in the sweep\n", label);
        return;
    }
    bisect_n_sub_switch(summary->worst_switch_klf_low, summary->worst_switch_klf_high,
                        lateral_analytic, theta_in, n_steps,
                        rain_mm_per_h, pet_mm_per_h, &klf_minus, &klf_plus);
    evaluate_tangent_nsub(klf_minus, theta_in, n_steps, 0, lateral_analytic,
                          rain_mm_per_h, pet_mm_per_h, y_minus, dy_minus, n_sub_minus);
    evaluate_tangent_nsub(klf_plus, theta_in, n_steps, 0, lateral_analytic,
                          rain_mm_per_h, pet_mm_per_h, y_plus, dy_plus, n_sub_plus);

    printf("%s: worst switch at k_lf = %.10e m/h\n", label, klf_minus);
    for (int i_step = 0; i_step < n_steps; i_step++) {
        if (n_sub_minus[i_step] != n_sub_plus[i_step]) {
            printf("   hour %d: n_sub %d -> %d   (rain %.3f mm/h)\n",
                   i_step + 1, n_sub_minus[i_step], n_sub_plus[i_step],
                   rain_mm_per_h[i_step]);
        }
    }
    printf("   %-22s %14s %14s %14s %12s\n", "output", "jump J", "AD left",
           "AD right", "J/(AD*k)");
    {
        double jump = total_lateral_from_vector(y_plus) - total_lateral_from_vector(y_minus);
        double ad_left = total_lateral_from_vector(dy_minus);
        double ad_right = total_lateral_from_vector(dy_plus);
        printf("   %-22s %14.6e %14.7e %14.7e %12.3e\n", "lateral_total_m",
               jump, ad_left, ad_right, jump / (ad_left * klf_minus));
    }
    for (int j = 0; j < 5; j++) {
        int k = show_index[j];
        double jump = y_plus[k] - y_minus[k];
        double equivalent = 0.0;

        if (fabs(dy_minus[k]) > TINY_DERIVATIVE_SCALE) {
            equivalent = jump / (dy_minus[k] * klf_minus);
        }
        printf("   %-22s %14.6e %14.7e %14.7e %12.3e\n", output_label[k],
               jump, dy_minus[k], dy_plus[k], equivalent);
    }
}

static void experiment_6_exponential_lateral(const char *output_dir,
                                             const char *forcing_path)
{
    static double rain_mm_per_h[E4_N_STEPS_H];
    static double pet_mm_per_h[E4_N_STEPS_H];
    static int n_sub_scratch[MAX_STEPS];
    const double single_step_klf[3] = {1.0e-3, 5.0e-2, 3.0e-1};
    const double compare_klf[3] = {BASELINE_KLF_M_PER_H, 1.0e-3, 1.0e-2};
    const char *scheme_label[2] = {"Euler", "exponential"};
    double theta_in[NDISC];
    double theta_fc;
    double zero_forcing[1];
    int n_steps;
    SweepSummary adaptive_summary;
    SweepSummary fixed_summary;
    FILE *fp;

    theta_fc = klf_experiment_theta_fc();
    baseline_entering_state(theta_in);
    zero_forcing[0] = 0.0;

    printf("\n====================================================================\n");
    printf("E6  EXACT EXPONENTIAL LATERAL REMOVAL versus FORWARD EULER\n");
    printf("====================================================================\n");

    /* --- E6a one timestep, baseline state, no forcing --- */
    printf("E6a  one timestep, baseline state, no rain, no PET (n_sub = 1)\n");
    printf("     analytic exponential check: d(lat_i)/dk = (lat_i/k) * a_i/(exp(a_i)-1),\n");
    printf("     a_i = k * dt / (dz_i * (theta_sat - theta_fc))\n");
    for (int i_klf = 0; i_klf < 3; i_klf++) {
        double klf = single_step_klf[i_klf];
        double h = klf * 1.0e-6;
        double y_euler[N_OUTPUTS];
        double dy_euler[N_OUTPUTS];
        double y_exp[N_OUTPUTS];
        double dy_exp[N_OUTPUTS];
        double y_plus[N_OUTPUTS];
        double y_minus[N_OUTPUTS];

        evaluate_tangent_nsub(klf, theta_in, 1, 0, 0, zero_forcing, zero_forcing,
                              y_euler, dy_euler, n_sub_scratch);
        evaluate_tangent_nsub(klf, theta_in, 1, 0, 1, zero_forcing, zero_forcing,
                              y_exp, dy_exp, n_sub_scratch);
        evaluate_primal_nsub(klf + h, theta_in, 1, 0, 1, zero_forcing, zero_forcing,
                             y_plus, n_sub_scratch);
        evaluate_primal_nsub(klf - h, theta_in, 1, 0, 1, zero_forcing, zero_forcing,
                             y_minus, n_sub_scratch);

        printf("\n   k_lf = %.3e m/h   n_sub = %d\n", klf, n_sub_scratch[0]);
        printf("   %-6s %14s %14s %14s %14s %14s %10s\n", "disc",
               "Euler lat", "Euler AD", "exp lat", "exp AD", "exp analytic", "FD-AD rel");
        for (int i_disc = 0; i_disc < NDISC; i_disc++) {
            double a = klf * 1.0 / (disc_thickness_m[i_disc] *
                                    (THETA_SAT_M3_PER_M3 - theta_fc));
            double analytic = (y_exp[i_disc] / klf) * a / (exp(a) - 1.0);
            double fd = (y_plus[i_disc] - y_minus[i_disc]) / (2.0 * h);

            printf("   %-6d %14.7e %14.7e %14.7e %14.7e %14.7e %10.2e\n", i_disc,
                   y_euler[i_disc], dy_euler[i_disc], y_exp[i_disc], dy_exp[i_disc],
                   analytic, relative_difference_from_ad(fd, dy_exp[i_disc]));
        }
        printf("   derivative volume residual: Euler %.2e   exponential %.2e\n",
               derivative_volume_residual(dy_euler), derivative_volume_residual(dy_exp));
    }

    n_steps = read_forcing_csv(forcing_path, E4_N_STEPS_H, rain_mm_per_h, pet_mm_per_h);
    if (n_steps <= 0) {
        printf("\nE6b-d skipped: could not read forcing file %s\n", forcing_path);
        return;
    }

    /* --- E6b exponential, adaptive n_sub sweep --- */
    printf("\nE6b-c  same %d-hour observed forcing and k_lf sweep as E5\n", n_steps);
    fp = open_csv(output_dir, "e6_sweep_exponential_adaptive.csv");
    sweep_klf(0, 1, theta_in, n_steps, rain_mm_per_h, pet_mm_per_h, fp, &adaptive_summary);
    fclose(fp);
    print_sweep_summary("E6b EXPONENTIAL lateral, ADAPTIVE n_sub", &adaptive_summary);
    print_switch_jumps("E6b", 1, &adaptive_summary, theta_in, n_steps,
                       rain_mm_per_h, pet_mm_per_h);

    /* --- E6c exponential, fixed n_sub sweep --- */
    fp = open_csv(output_dir, "e6_sweep_exponential_fixed12.csv");
    sweep_klf(E5_FIXED_N_SUB, 1, theta_in, n_steps, rain_mm_per_h, pet_mm_per_h,
              fp, &fixed_summary);
    fclose(fp);
    print_sweep_summary("E6c EXPONENTIAL lateral, FIXED n_sub = 12", &fixed_summary);

    /* --- E6d sensitivity of values and derivatives to n_sub, both schemes --- */
    printf("\nE6d  adaptive versus fixed n_sub = 12, for each lateral scheme\n");
    printf("     (rel = (fixed12 - adaptive)/|adaptive|; smaller means less\n");
    printf("      dependence on the time discretization)\n");
    printf("%-9s %-18s %-12s %15s %15s %10s %10s\n", "k_lf", "output", "scheme",
           "AD adaptive", "AD fixed12", "AD rel", "value rel");
    for (int i_klf = 0; i_klf < 3; i_klf++) {
        double klf = compare_klf[i_klf];

        for (int scheme = 0; scheme < 2; scheme++) {
            double y_a[N_OUTPUTS];
            double dy_a[N_OUTPUTS];
            double y_f[N_OUTPUTS];
            double dy_f[N_OUTPUTS];
            double key_value_a[3];
            double key_value_f[3];
            double key_ad_a[3];
            double key_ad_f[3];
            const char *key_name[3] = {"lateral_total_m", "percolation_m", "theta_out[0]"};

            evaluate_tangent_nsub(klf, theta_in, n_steps, 0, scheme,
                                  rain_mm_per_h, pet_mm_per_h, y_a, dy_a, n_sub_scratch);
            evaluate_tangent_nsub(klf, theta_in, n_steps, E5_FIXED_N_SUB, scheme,
                                  rain_mm_per_h, pet_mm_per_h, y_f, dy_f, n_sub_scratch);

            key_value_a[0] = total_lateral_from_vector(y_a);
            key_value_f[0] = total_lateral_from_vector(y_f);
            key_ad_a[0] = total_lateral_from_vector(dy_a);
            key_ad_f[0] = total_lateral_from_vector(dy_f);
            key_value_a[1] = y_a[IDX_PERC];
            key_value_f[1] = y_f[IDX_PERC];
            key_ad_a[1] = dy_a[IDX_PERC];
            key_ad_f[1] = dy_f[IDX_PERC];
            key_value_a[2] = y_a[IDX_THETA_0];
            key_value_f[2] = y_f[IDX_THETA_0];
            key_ad_a[2] = dy_a[IDX_THETA_0];
            key_ad_f[2] = dy_f[IDX_THETA_0];

            for (int j = 0; j < 3; j++) {
                printf("%-9.1e %-18s %-12s %15.8e %15.8e %10.2e %10.2e\n",
                       klf, key_name[j], scheme_label[scheme],
                       key_ad_a[j], key_ad_f[j],
                       relative_difference_from_ad(key_ad_f[j], key_ad_a[j]),
                       relative_difference_from_ad(key_value_f[j], key_value_a[j]));
            }
            printf("%-9.1e %-18s %-12s derivative volume residual: adaptive %.2e  fixed12 %.2e\n",
                   klf, "", scheme_label[scheme],
                   derivative_volume_residual(dy_a), derivative_volume_residual(dy_f));
        }
    }
}


/*
 * Consistency check run before any experiment: the primal wrapper (linked
 * against libsoil.a) and the Tapenade tangent (which carries its own copy
 * of the kernel) must produce the same primal values.  A mismatch means
 * the library and the generated code were built from different sources
 * or struct layouts, and every finite-difference result would be wrong.
 */
static int primal_and_tangent_values_agree(void)
{
    const double klf_check_m_per_h = 1.0e-3;
    double theta_in[NDISC];
    double rain_mm_per_h[1];
    double pet_mm_per_h[1];
    double y_primal[N_OUTPUTS];
    double y_tangent[N_OUTPUTS];
    double dy_tangent[N_OUTPUTS];
    int n_sub_primal[1];
    int n_sub_tangent[1];
    int all_agree;

    baseline_entering_state(theta_in);
    rain_mm_per_h[0] = 2.0;
    pet_mm_per_h[0] = 0.3;

    evaluate_primal(klf_check_m_per_h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                    y_primal, n_sub_primal);
    evaluate_tangent(klf_check_m_per_h, theta_in, 1, rain_mm_per_h, pet_mm_per_h,
                     y_tangent, dy_tangent, n_sub_tangent);

    all_agree = 1;
    if (n_sub_primal[0] != n_sub_tangent[0]) all_agree = 0;
    for (int k = 0; k < N_OUTPUTS; k++) {
        double scale = fabs(y_tangent[k]);

        if (scale < 1.0e-12) scale = 1.0e-12;
        if (fabs(y_primal[k] - y_tangent[k]) > 1.0e-12 * scale) {
            printf("CONSISTENCY FAILURE %-22s primal %.17e  tangent %.17e\n",
                   output_label[k], y_primal[k], y_tangent[k]);
            all_agree = 0;
        }
    }
    if (!all_agree) {
        printf("ERROR: primal (libsoil.a) and Tapenade tangent disagree on primal values.\n");
        printf("       libsoil.a is probably stale.  Rebuild with: make veryclean ; make\n");
    }
    return all_agree;
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    const char *output_dir = ".";
    const char *forcing_path = "forcing/rain_pet_example.csv";

    const char *experiments = "123456";

    if (argc > 1) output_dir = argv[1];
    if (argc > 2) forcing_path = argv[2];
    if (argc > 3) experiments = argv[3];

    if (!primal_and_tangent_values_agree()) return 2;
    printf("Consistency check passed: primal and tangent primal values agree.\n");

    if (strchr(experiments, '1') != NULL) experiment_1_single_step(output_dir);
    if (strchr(experiments, '2') != NULL) experiment_2_storage_cap_kink(output_dir);
    if (strchr(experiments, '3') != NULL) experiment_3_drydown(output_dir);
    if (strchr(experiments, '4') != NULL) experiment_4_observed_forcing(output_dir, forcing_path);
    if (strchr(experiments, '5') != NULL) experiment_5_n_sub_boundaries(output_dir, forcing_path);
    if (strchr(experiments, '6') != NULL) experiment_6_exponential_lateral(output_dir, forcing_path);

    return 0;
}
