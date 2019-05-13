#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "graph.h"
#include "trigonometry.h"
#include "io.h"
#include "opt.h"
#include "core_l.h"
#include "core_b.h"
#include "core_vr.h"
#include "math.h"
#include "utils.h"
#include "unicore.h"
#include "mem.h"


#define OUTPUT_RESULT_FILENAME  \
        "result.txt"

#define OUTPUT_UNFRESULT_FILENAME  \
        "unfresult.txt"

#define OUTPUT_RESULT_LINE      \
        "+---------------------------------------------+"

#define PRINT_OUTPUT_LINE(f)    \
        fprintf((f), "%s\n", OUTPUT_RESULT_LINE)


#define CHECK_FILE_AND_RET(fdesc, name)                 \
        do {                                            \
                if ((fdesc) == NULL) {                  \
                        PRINT_IO_OPEN_ERROR(name);      \
                        return;                         \
                }                                       \
        } while (0)

parameter_t g_ptable[] = {
        { "u" }
      , { "v" }
      , { "w" }
//      , { "o_0" }
      , { "A" }
        // MUST BE LAST:
      , { NULL }
};

#define MAX_G_GRAPH_LEN         8
static char __g_graph_buffer[MAX_G_GRAPH_LEN];

static char *__get_name_by_idx(const unsigned int idx)
{
        const unsigned int known = BETA_QTY + 1;
        if (idx < known) return g_ptable[idx].name;

        sprintf(__g_graph_buffer, "th[%u]", idx - known + 2);
        return __g_graph_buffer;
}

static void dump_unfriendly_result(const opt_t *opt)
{
        FILE *fout = fopen(OUTPUT_UNFRESULT_FILENAME, "w");
        CHECK_FILE_AND_RET(fout, OUTPUT_UNFRESULT_FILENAME);

        fprintf(fout, "%0.3lf %0.3lf\n",
                        opt->r_0,
                        opt->dr_0);
        //fprintf(fout, "%0.3lf\n",
        //                sqrt(opt->sq / (opt->size + opt->s.size + 1)));

        //fprintf(fout, "%d\n", opt->s.size - BETA_QTY);
        unsigned int i;
        for (i = 0; i < opt->s.size; ++i) {
                fprintf(fout, "%0.3f %0.3f %0.2lf\n",
                        opt->s.data[i],
                        opt->bounds[i].l,
                        fabs(opt->bounds[i].l / opt->s.data[i]));
        }
        fclose(fout);
}

void dump_vr_solution(const opt_t *opt)
{
        FILE *fout = fopen(OUTPUT_RESULT_FILENAME, "w");
        CHECK_FILE_AND_RET(fout, OUTPUT_RESULT_FILENAME);

        PRINT_OUTPUT_LINE(fout);
        fprintf(fout, "Result for APOGEE-RC dataset by %u obj:\n",
                        opt->size);
        PRINT_OUTPUT_LINE(fout);
        fprintf(fout, "R_0: \t%0.3lf \t+%0.3lf\n",
                        opt->r_0, opt->dr_0);
        fprintf(fout, "SD: \t%0.3lf\n",
                        opt->sq);
        PRINT_OUTPUT_LINE(fout);
        unsigned int i;
        for (i = 0; i < opt->s.size; ++i) {
                fprintf(fout, "%s:\t%6.3f\t(pm %0.3f)\t%0.2lf\n",
                        __get_name_by_idx(i),
                        opt->s.data[i],
                        opt->bounds[i].l,
                        fabs(opt->bounds[i].l / opt->s.data[i]));
        }
        PRINT_OUTPUT_LINE(fout);
        fclose(fout);


}

void dump_result(const opt_t *opt)
{
        FILE *fout = fopen(OUTPUT_RESULT_FILENAME, "w");
        CHECK_FILE_AND_RET(fout, OUTPUT_RESULT_FILENAME);

        PRINT_OUTPUT_LINE(fout);
        fprintf(fout, "Result for APOGEE-RC dataset by %u obj:\n",
                        opt->size);
        PRINT_OUTPUT_LINE(fout);
        fprintf(fout, "R_0: \t%0.3lf \t+%0.3lf\n",
                        opt->r_0, opt->dr_0);
        fprintf(fout, "\t\t-%0.3lf\n", opt->dr_0);
        fprintf(fout, "SD: \t%0.3lf\n",
                        sqrt(opt->sq / (3 * opt->size - opt->s.size - 1)));
        PRINT_OUTPUT_LINE(fout);
        unsigned int i;
#ifdef DEBUG
        printf("%s: opt->s.size = %u\n",
                        __func__, opt->s.size);
        for (i = 0; i < g_ptable[i].name != NULL; ++i) {
                printf("%s: [%u] = %lf\n", __func__,
                                i, opt->s.data[i]);
        }
#endif
        for (i = 0; i < opt->s.size; ++i) {
                fprintf(fout, "%s:\t%6.3f\t(pm %0.3f)\t%0.2lf\n",
                        __get_name_by_idx(i),
                        opt->s.data[i],
                        opt->bounds[i].l,
                        fabs(opt->bounds[i].l / opt->s.data[i]));
        }
        PRINT_OUTPUT_LINE(fout);
        fclose(fout);

        dump_unfriendly_result(opt);
}


double get_point_by_uni_solution(const opt_t *solution,
                                 const double r)
{
        double theta = r * solution->s.data[BETA_QTY] -
                                        2 * solution->s.data[BETA_QTY + 1] * (r - solution->r_0);
        unsigned int i;
        for (i = BETA_QTY + 2; i < solution->s.size; ++i) {
                theta += solution->s.data[i] / dv_factorial(i - BETA_QTY) *
                                        pow_double(r - solution->r_0, i - BETA_QTY);
        }

        return theta;
}

double get_point_by_vr_solution(const opt_t *solution,
                                const double r)
{
        double theta = r * (OMEGA_SUN - solution->s.data[V_P] / solution->r_0) -
                                        2 * solution->s.data[BETA_QTY] * (r - solution->r_0);
        unsigned int i;
        for (i = BETA_QTY + 1; i < solution->s.size; ++i) {
                theta += solution->s.data[i] / dv_factorial(i - BETA_QTY + 1) *
                                        pow_double(r - solution->r_0, i - BETA_QTY + 1);
        }

        return theta;
}

double get_point_by_b_solution(const opt_t *solution,
                               const double r)
{
        double theta = r * (OMEGA_SUN - solution->s.data[V_P] / solution->r_0) +
                                        2 * solution->s.data[BETA_QTY] * (r - solution->r_0);
        unsigned int i;
        for (i = BETA_QTY + 1; i < solution->s.size; ++i) {
                theta -= solution->s.data[i] / dv_factorial(i - BETA_QTY + 1) *
                                        pow_double(r - solution->r_0, i - BETA_QTY + 1);
        }

        return theta;
}

double get_point_by_l_solution(const opt_t *solution,
                               const double r)
{
        double theta = r * solution->s.data[W_P] -
                                        2 * solution->s.data[BETA_QTY] * (r - solution->r_0);
        unsigned int i;
        for (i = BETA_QTY + 1; i < solution->s.size; ++i) {
                theta += solution->s.data[i] / dv_factorial(i - BETA_QTY + 1) *
                                        pow_double(r - solution->r_0, i - BETA_QTY + 1);
        }

        return theta;
}



double get_c_point_by_part_solution(const opt_t *solution,
                                    const double r,
                                    unsigned int type,
                                    const double omega_0,
                                    unsigned int A_idx)
{
        double theta = 2 * solution->s.data[A_idx] * (r - solution->r_0);
        double theta_s = 0;

        unsigned int i;
        for (i = A_idx + 1; i < solution->s.size; ++i) {
                theta_s += solution->s.data[i] / dv_factorial(i - A_idx + 1) * 
                                        pow_double(r - solution->r_0, i - A_idx + 1);
        }

        double theta_omega = r * omega_0;

        return type != B_PART ? -theta + theta_s + theta_omega
                              : theta - theta_s + theta_omega;
}




static double get_mu_sun(const apogee_rc_t *line,
                         const linear_eq_solve_t *solution,
                         double (*f_beta)(const apogee_rc_t *,
                                          beta_ord_t),
                        unsigned int A_idx)
{
        unsigned int i;
        double res = 0;

        for (i = 0; i < A_idx; ++i) {
                res += f_beta(line, i) * solution->data[i];
        }

        return res;
}



/**
 * Object to Theta-R diagramm
 */
static double obs_theta_R_by_l(const opt_t *solution,
                               const apogee_rc_t *line,
                               const double omega_0)
{
        double R = get_R_distance(line, solution->r_0);
        double mu_l = get_mu_sun(line, &solution->s, core_l_get_beta_n, BETA_QTY);
        return ((line->pm_l - mu_l) /
                        (solution->r_0 * line->cos_l / line->dist - line->cos_b) + omega_0) * R;
}


static double obs_theta_R_by_b(const opt_t *solution,
                               const apogee_rc_t *line,
                               const double omega_0)
{
        double R = get_R_distance(line, solution->r_0);
        double mu_b = get_mu_sun(line, &solution->s, core_b_get_beta_n, BETA_QTY);
        return ((-line->pm_b + mu_b) / (solution->r_0 * line->sin_l * line->sin_b) * line->dist + omega_0) * R;
}


static double obs_theta_R_by_vr(const opt_t *solution,
                                const apogee_rc_t *line,
                                const double omega_0)
{
        double R = get_R_distance(line, solution->r_0);
        double vr_sun = get_mu_sun(line, &solution->s, __core_vr_get_beta_n, BETA_QTY_FIX);
        return ((line->v_helio - vr_sun - __core_vr_get_beta_n(line, THIRD)) / (solution->r_0 * line->sin_l * line->cos_b) + omega_0) * R;
}

static double obs_theta_R_by_vr_free(const opt_t *solution,
                                     const apogee_rc_t *line,
                                     const double omega_0)
{
        double R = get_R_distance(line, solution->r_0);
        double vr_sun = get_mu_sun(line, &solution->s, core_vr_get_beta_n, BETA_QTY);
        return ((line->v_helio - vr_sun) / (solution->r_0 * line->sin_l * line->cos_b) + omega_0) * R;
}

double theta_by_R_vr(const opt_t *solution, const double r)
{
        return 0;
}

/**
 * Generic dumper for objects.
 *
 * @type: see to enum VR_PART, B_PART, L_PART
 */
void dump_objects_theta_R(const apogee_rc_table_t *table,
                          opt_t *solution,
                          unsigned int type,
                          const char *filename)
{
        unsigned int size = table->size;
        unsigned int i;

        FILE *fout = fopen(filename, "w");
        CHECK_FILE_AND_RET(fout, filename);

        double (*f)(const opt_t *, const apogee_rc_t *, const double);

        if (type == VR_PART) {
                f = obs_theta_R_by_vr;
        }

        if (type == L_PART) {
                //solution->s.data[BETA_QTY - 1] = 0;
                f = obs_theta_R_by_l;
        }

        if (type == B_PART) {
                f = obs_theta_R_by_b;
        }

        if (type == TOTAL_QTY) {
                f = obs_theta_R_by_vr_free;
        }

        for (i = 0; i < size; ++i) {
                double R = get_R_distance(&table->data[i], solution->r_0);
                double theta = f(solution, &table->data[i], table->omega_0);
                fprintf(fout, "%lf %lf\n", R, theta);
        }

        fclose(fout);
}

void dump_part_rotation_curve(const opt_t *solution,
                              unsigned int type,
                              const char *filename,
                              const double omega_0)
{
        double R = ROTC_LOWER_BOUND;
        FILE *fout = fopen(filename, "w");
        CHECK_FILE_AND_RET(fout, filename);
        unsigned int A_idx = type == VR_MODE ? BETA_QTY_FIX
                                             : BETA_QTY;

        while (R < ROTC_UPPER_BOUND) {
                double theta = get_c_point_by_part_solution(solution,
                                                            R, type,
                                                            omega_0,
                                                            A_idx);
                fprintf(fout, "%lf %lf\n", R, theta);
                R += ROTC_STEP_R;
        }

        fclose(fout);
}

void dump_uni_rotation_curve(const rot_curve_t *curve, const unsigned int size)
{
        unsigned int i;

        FILE *fout = fopen(RC_OUT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(fout, RC_OUT_FILE_NAME);

        for (i = 0; i < size; ++i) {
                fprintf(fout, "%lf %lf %lf %lf\n",
                                curve[i].r, curve[i].theta, curve[i].theta_max, curve[i].theta_min);
        }
        fclose(fout);
}

void dump_rotation_curve_vr(const opt_t *solution)
{
        FILE *fout = fopen(RC_OUT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(fout, RC_OUT_FILE_NAME);

        FILE *sout = fopen(SUN_POINT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(sout, SUN_POINT_FILE_NAME);

#ifdef DEBUG
        printf("%s: enrty\n", __func__);
#endif
        unsigned int i;
        double r;

        r = ROTC_LOWER_BOUND;
        while (r < ROTC_UPPER_BOUND) {
                double theta = r * (OMEGA_SUN - solution->s.data[V_P] / solution->r_0) -
                                        2 * solution->s.data[A_P] * (r - solution->r_0);
                for (i = BETA_QTY + 1; i < solution->s.size; ++i) {
                        theta += solution->s.data[i] / dv_factorial(i - BETA_QTY + 1) * 
                                        pow_double(r - solution->r_0, i - BETA_QTY + 1);
                }

                fprintf(fout, "%lf %lf\n", r, theta);
                r += ROTC_STEP_R;
        }

        fprintf(sout, "%lf %lf\n", solution->r_0, solution->r_0 * OMEGA_SUN);

        fclose(sout);
        fclose(fout);
}

void dump_rotation_curve_b(const opt_t *solution)
{
        FILE *fout = fopen(RC_OUT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(fout, RC_OUT_FILE_NAME);

        FILE *sout = fopen(SUN_POINT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(sout, SUN_POINT_FILE_NAME);

#ifdef DEBUG
        printf("%s: enrty\n", __func__);
#endif
        unsigned int i;
        double r;

        r = ROTC_LOWER_BOUND;
        while (r < ROTC_UPPER_BOUND) {
                double theta = r * (OMEGA_SUN - solution->s.data[V_P] / solution->r_0) +
                                        2 * solution->s.data[A_P] * (r - solution->r_0);
                for (i = BETA_QTY + 1; i < solution->s.size; ++i) {
                        theta -= solution->s.data[i] / dv_factorial(i - BETA_QTY + 1) *
                                        pow_double(r - solution->r_0, i - BETA_QTY + 1);
                }

                fprintf(fout, "%lf %lf\n", r, theta);
                r += ROTC_STEP_R;
        }

        fprintf(sout, "%lf %lf\n", solution->r_0, solution->r_0 * OMEGA_SUN);

        fclose(sout);
        fclose(fout);
}

void dump_rotation_curve_l(const opt_t *solution)
{
        FILE *fout = fopen(RC_OUT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(fout, RC_OUT_FILE_NAME);

        FILE *sout = fopen(SUN_POINT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(sout, SUN_POINT_FILE_NAME);

#ifdef DEBUG
        printf("%s: enrty\n", __func__);
#endif
        unsigned int i;
        double r;

        r = ROTC_LOWER_BOUND;
        while (r < ROTC_UPPER_BOUND) {
                double theta = r * solution->s.data[W_P] -
                                        2 * solution->s.data[A_P] * (r - solution->r_0);
                for (i = BETA_QTY + 1; i < solution->s.size; ++i) {
                        theta += solution->s.data[i] / dv_factorial(i - BETA_QTY + 1) *
                                        pow_double(r - solution->r_0, i - BETA_QTY + 1);
                }

                fprintf(fout, "%lf %lf\n", r, theta);
                r += ROTC_STEP_R;
        }

        fprintf(sout, "%lf %lf\n", solution->r_0, solution->r_0 * (solution->s.data[W_P] + solution->s.data[V_P] / solution->r_0));

        fclose(sout);
        fclose(fout);
}


static double *get_sorted_r(const iteration_storage_t *sorted_st, const size_t size)
{
        /**
	 * BUG: double *sorted_r = (double *)dv_dalloc(sizeof(double), size);
	 *
	 * Not valid code: causes crash by access memory failure.
	 */
        double *sorted_r = (double *)calloc(size, sizeof(double));
        unsigned int i = 0;
        for (i = 0; i < size; ++i) {
		sorted_r[i] = sorted_st[i].r;
        }

        return sorted_r;
}

void dump_averages(iteration_storage_t *st, opt_t *solution, averages_mode_t mode)
{
        average_res_t *a;

        sort_iteration_storage_by_r(st, solution->size);
        double *sorted_r = get_sorted_r(st, solution->size);

        FILE *aout = fopen(AVERAGE_R_FILE_NAME, "w");
        CHECK_FILE_AND_RET(aout, AVERAGE_R_FILE_NAME);

        double r = st[0].r;
        int estimate_counter = (int)solution->size;
        int left_bound = 0;
        int i_counter = 1;

        while (estimate_counter > 0) {
                /* Setting for size */
                int size = R_INTO_MIDDLE(r) ? AVERAGE_COUNT_BASE
                                            : AVERAGE_COUNT_EDGE;

                if (IS_LAST_STEP(estimate_counter, size)) {
                        size = estimate_counter;
                }

                /* Get & print */
                a = get_average_theta(&st[left_bound], solution, (unsigned int)size);
                double median_r = get_median(&sorted_r[left_bound], (unsigned int)size);
                fprintf(aout, "%d %d %lf %lf %lf %lf %lf\n",
                                ++i_counter, size,
                                median_r,
                                fabs(median_r - sorted_r[left_bound]),
                                fabs(median_r - sorted_r[left_bound + size - 1]),
                                a->theta, a->err
                                );


                /*  Next step prepare */
                r = sorted_r[left_bound + size - 1];
                left_bound += size;
                estimate_counter -= size;
        }

        fclose(aout);
}

void dump_R0_theta_ellips(const opt_t **res,
                          unsigned int size,
                          const opt_t *solution)
{
        unsigned int i;

        FILE *aout = fopen(R0THETA0_FILE_NAME, "w");
        CHECK_FILE_AND_RET(aout, R0THETA0_FILE_NAME);

        for (i = 0; i < size; ++i) {
                fprintf(aout, "%lf %lf\n",
                        res[i]->r_0, res[i]->s.data[BETA_QTY] * res[i]->r_0);
        }

        FILE *cout = fopen(R0THETA0_MAIN_FILE_NAME, "w");
        CHECK_FILE_AND_RET(aout, R0THETA0_MAIN_FILE_NAME);

        fprintf(cout, "%lf %lf", solution->r_0, 
                                 solution->r_0 * solution->s.data[BETA_QTY]);
        fclose(cout);
        fclose(aout);
}

void dump_core_l_solution(const opt_t *solution)
{
        printf("L Core Solution:\n");
        printf("u_sun: %lf pm %lf\n", solution->s.data[0],
                                      solution->bounds[0].l);
        printf("v_sun: %lf pm %lf\n", solution->s.data[1],
                                      solution->bounds[1].l);
        printf("omega_0: %lf pm %lf\n", solution->s.data[2],
                                      solution->bounds[2].l);
        printf("A: %lf pm %lf\n", solution->s.data[3],
                                      solution->bounds[3].l);
        printf("pm_l SD: %lf\n", solution->sq);
        printf("size: %d\n", solution->size);
        printf("---------------\n");
}

void dump_core_vr_solution(const opt_t *solution)
{
        printf("Vr Core Solution:\n");
        printf("R_0: %lf \n", solution->r_0);
        printf("u_sun: %lf \n", solution->s.data[0]);
        printf("v_sun: %lf \n", solution->s.data[1]);
        printf("A: %lf \n", solution->s.data[2]);
        printf("SD: %lf\n", solution->sq);
        printf("size: %d\n", solution->size);
        printf("---------------\n");
}

void dump_united_solution(const opt_t *solution)
{
        printf("United Solution:\n");
        printf("R_0: %lf\n", solution->r_0);
        printf("u_sun: %lf pm %lf\n", solution->s.data[0],
                                      solution->bounds[0].l);
        printf("v_sun: %lf pm %lf\n", solution->s.data[1],
                                      solution->bounds[1].l);
        printf("w_sun: %lf pm %lf\n", solution->s.data[2],
                                      solution->bounds[2].l);
        printf("omega_0: %lf pm %lf\n", solution->s.data[3],
                                      solution->bounds[3].l);
        printf("A: %lf pm %lf\n", solution->s.data[4],
                                      solution->bounds[4].l);
        printf("khi_sq: %lf\n", solution->sq);
        printf("N_free: %d\n", 3 * solution->size - solution->s.size - 1);
        printf("---------------\n");
}

void dump_united_solution_points(const opt_t *solution)
{
        printf("United Solution:\n");
        printf("R_0: %lf\n", solution->r_0);
        printf("u_sun: %lf\n", solution->s.data[0]);
        printf("v_sun: %lf\n", solution->s.data[1]);
        printf("w_sun: %lf\n", solution->s.data[2]);
        printf("omega_0: %lf\n", solution->s.data[3]);
        printf("A: %lf\n", solution->s.data[4]);
        printf("khi_sq: %lf\n", solution->sq);
        printf("N_free: %d\n", 3 * solution->size - solution->s.size - 1);
        printf("---------------\n");
}

void dump_table_parameters(const apogee_rc_table_t *table,
                           const opt_t *solution)
{
        printf("SD[VR_PART] = %lf\n", sqrt(table->sigma[VR_PART]));
        printf("SD[L_PART] = %lf\n", sqrt(table->sigma[L_PART]));
        printf("SD[B_PART] = %lf\n", sqrt(table->sigma[B_PART]));
}

void dump_core_b_solution(const opt_t *solution)
{
        printf("B Core Solution:\n");
        printf("u_sun: %lf pm %lf\n", solution->s.data[0], 
                                      solution->bounds[0].l);
        printf("v_sun: %lf pm %lf\n", solution->s.data[1],
                                      solution->bounds[1].l);
        printf("w_sun: %lf pm %lf\n", solution->s.data[2],
                                      solution->bounds[2].l);
        printf("A: %lf pm %lf\n", solution->s.data[3],
                                      solution->bounds[3].l);
        printf("pm_b SD: %lf\n", solution->sq);
        printf("size: %d\n", solution->size);
        printf("---------------\n");
}

void dump_background(const iteration_storage_t *st,
                     const opt_t *solution,
                     const int b_count)
{
        FILE *dout = fopen("background.txt", "w");
        CHECK_FILE_AND_RET(dout, "background.txt");

        average_res_t *a;
        double *sorted_r = get_sorted_r(st, solution->size);

        double r = st[0].r;
        int estimate_counter = (int)solution->size;
        int left_bound = 0;
        int size = (int)solution->size / b_count;

        while (estimate_counter > 0) {
                /* Setting for size */

                if (IS_LAST_STEP(estimate_counter, size)) {
                        size = estimate_counter;
                }

                /* Get & print */
                a = get_average_theta(&st[left_bound], solution, (unsigned int)size);

                fprintf(dout, "%lf %lf %lf %lf %lf %d %lf\n",
                               get_median(&sorted_r[left_bound], (unsigned int)size),
                               sorted_r[left_bound],
                               sorted_r[left_bound + size - 1],
                               a->sd,
                               a->err,
                               size,
                               sqrt(solution->sq / (solution->size + solution->s.size + 1)));

                /*  Next step prepare */
                r = sorted_r[left_bound + size - 1];
                left_bound += size;
                estimate_counter -= size;
        }

#if 0 // dummy now:
        FILE *fout = fopen("bk_sd.txt", "w");
        CHECK_FILE_AND_RET(fout, "bk_sd.txt");
        fprintf(fout, "%lf\n", 0);

        fclose(fout);
#endif
        fclose(dout);
}



void dump_rand_test(const double *array,
                    const dsize_t size)
{
        FILE *fout = fopen("rand.txt", "w");
        CHECK_FILE_AND_RET(fout, "rand.txt");

        dsize_t i;
        for (i = 0; i < size; i += 2) {
                fprintf(fout, "%0.7lf \t%0.7lf\n",
                                array[i], array[i + 1]);

        }

        fclose(fout);
}

static double cos_beta_uni(const apogee_rc_t *line,
                           const double r_0)
{
        return (r_0 - line->dist * cos(line->b) * cos(line->l)) / get_R_distance(line, r_0);
}

static double sin_beta_uni(const apogee_rc_t *line,
                           const double r_0)
{
        return line->dist * cos_beta_uni(line, r_0) * sin(line->l) / get_R_distance(line, r_0);
}

static double get_v_l(const apogee_rc_t *line)
{
        return K_PM * line->dist * line->pm_l * cos(line->b);
}

static double get_v_b(const apogee_rc_t *line)
{
        return K_PM * line->dist * line->pm_b;
}

static double get_u_g(const apogee_rc_t *line,
                      const opt_t *solution)
{
        return solution->s.data[0] + cos(line->l) * (line->v_helio * cos(line->b) - 
                        get_v_b(line) * sin(line->b)) - get_v_l(line) * sin(line->l);
}

static double get_v_g(const apogee_rc_t *line,
                      const opt_t *solution)
{
        return solution->s.data[BETA_QTY] * solution->r_0 + solution->s.data[1] + sin(line->l) * (line->v_helio * cos(line->b) - 
                        get_v_b(line) * sin(line->b)) + get_v_l(line) * cos(line->l);
}

static double get_theta_uni_obj(const apogee_rc_t *line,
                                const opt_t *solution)
{
        return get_v_g(line, solution) * cos_beta_uni(line, solution->r_0) +
                get_u_g(line, solution) * sin_beta_uni(line, solution->r_0);
}

void dump_uni_rotation_objs_named(const apogee_rc_table_t *table,
                                  const opt_t *solution,
                                  const char *filename)
{
        FILE *oout = fopen(filename, "w");
        CHECK_FILE_AND_RET(oout, filename);

        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                fprintf(oout, "%lf %lf\n",
                                get_R_distance(&table->data[i], solution->r_0),
                                get_theta_uni_obj(&table->data[i], solution));
        }
        fclose(oout);
}

void dump_uni_rotation_objs(const apogee_rc_table_t *table,
                            const opt_t *solution)
{
        dump_uni_rotation_objs_named(table,
                                     solution,
                                     DF_OUT_FILE_NAME);

        FILE *sout = fopen(SUN_POINT_FILE_NAME, "w");
        CHECK_FILE_AND_RET(sout, SUN_POINT_FILE_NAME);

        fprintf(sout, "%lf %lf\n", solution->r_0, solution->r_0 * solution->s.data[BETA_QTY] + solution->s.data[1]);

        fclose(sout);
}

void dump_line_xyz(const apogee_rc_t *line)
{
        // TODO: implement
        return;
}


void dump_objects_xyz_is(const iteration_storage_t *storage)
{
        // TODO: implement
        return;
}


void dump_objects_xyz(const apogee_rc_table_t *table,
                      const dsize_t size,
                      const char *filename)
{
        FILE *fout = fopen(filename, "w");
        CHECK_FILE_AND_RET(fout, filename);

        assert(table->size >= size);
        dsize_t i;
        for (i = 0; i < size; ++i) {
                point_t *p = get_point(&table->data[i]);
                fprintf(fout, "%0.7lf %0.7lf %0.7lf\n",
                                p->x, p->y, p->z);
        }

        fclose(fout);
}

void dump_table(const apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        const char *fname = cfg->dump_file_name;
        if (!strcmp(fname, cfg->input_file_name)) {
                printf("%s: [warning] input and dump files are equals!\n",
                                __func__);
        }

        FILE *fout = fopen(fname, "w");
        CHECK_FILE_AND_RET(fout, fname);

        dsize_t i;
        for (i = 0; i < table->size; ++i) {
                fprintf(fout, "%0.7lf %0.7lf %0.7lf %0.7lf %0.7lf %0.7lf\n",
                                rad_to_deg(table->data[i].l),
                                rad_to_deg(table->data[i].b),
                                table->data[i].v_helio,
                                table->data[i].dist,
                                table->data[i].dist,
                                table->data[i].dist);
        }

        fclose(fout);
}
