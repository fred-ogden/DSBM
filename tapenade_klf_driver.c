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
 * Usage:
 *   tapenade_klf_derivative_test  output_directory  forcing_csv_file
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

/* Run the primal wrapper and pack all dependents into one vector. */
static void evaluate_primal(double klf_m_per_h,
                            const double *theta_in,
                            int n_steps,
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

    dsbm_lateral_from_klf(klf_m_per_h, theta_in, n_steps,
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
static void evaluate_tangent(double klf_m_per_h,
                             const double *theta_in,
                             int n_steps,
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

    klf_tangent_from_tapenade(klf_m_per_h, theta_in, n_steps,
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

int main(int argc, char **argv)
{
    const char *output_dir = ".";
    const char *forcing_path = "forcing/rain_pet_example.csv";

    if (argc > 1) output_dir = argv[1];
    if (argc > 2) forcing_path = argv[2];

    experiment_1_single_step(output_dir);
    experiment_2_storage_cap_kink(output_dir);
    experiment_3_drydown(output_dir);
    experiment_4_observed_forcing(output_dir, forcing_path);

    return 0;
}
