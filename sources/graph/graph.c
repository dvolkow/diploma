#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "graph.h"
#include "io.h"
#include "opt.h"
#include "math.h"
#include "utils.h"


#define OUTPUT_RESULT_FILENAME  \
        "result.txt"

#define OUTPUT_RESULT_LINE      \
        "+---------------------------------------------+"

#define PRINT_OUTPUT_LINE(f)    \
        fprintf((f), "%s\n", OUTPUT_RESULT_LINE)


parameter_t g_ptable[] = {
        { "u" }
      , { "v" }
      , { "w" }
      , { "A" }
        // MUST BE LAST:
      , { NULL }
};

#define MAX_G_GRAPH_LEN         8
static char __g_graph_buffer[MAX_G_GRAPH_LEN];

static char *__get_name_by_idx(const int idx) 
{
        const int known = BETA_QTY + 1;
        if (idx < known) return g_ptable[idx].name;

        sprintf(__g_graph_buffer, "th[%d]", idx - known + 2);
        return __g_graph_buffer;
}

void dump_result(opt_t *opt, apogee_rc_table_t *table,
                  prec_t *p) 
{
        FILE *fout = fopen(OUTPUT_RESULT_FILENAME, "w"); 
        if (fout == NULL) {
                PRINT_IO_OPEN_ERROR(OUTPUT_RESULT_FILENAME);
                return;
        }

        PRINT_OUTPUT_LINE(fout);
        fprintf(fout, "Result for APOGEE-RC dataset by %lu obj:\n",
                        table->size);
        PRINT_OUTPUT_LINE(fout);
        fprintf(fout, "R_0: \t%lf \t+%lf\n",
                        opt->r_0, p->h - opt->r_0);
        fprintf(fout, "\t\t    \t-%lf\n", opt->r_0 - p->l);
        fprintf(fout, "SD: \t%lf\n",
                        sqrt(opt->sq / (table->size + opt->s.size + 1)));
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
                fprintf(fout, "%s: \t%f \t(pm %f) \n",
                        __get_name_by_idx(i),
                        opt->s.data[i],
                        opt->bounds[i].l);
        }
        fclose(fout);
}


/*
static double get_sun_shift(const apogee_rc_table_t *table, 
                                const opt_t *solution, const int idx)
{
}
*/

void dump_rotation_curve(iteration_storage_t *storage, opt_t *solution)
{
        FILE *fout = fopen(RC_OUT_FILE_NAME, "w");
        if (fout == NULL) {
                PRINT_IO_OPEN_ERROR(RC_OUT_FILE_NAME);
                return;
        }

        FILE *oout = fopen(DF_OUT_FILE_NAME, "w");
        if (oout == NULL) {
                PRINT_IO_OPEN_ERROR(RC_OUT_FILE_NAME);
                return;
        }

        FILE *sout = fopen(SUN_POINT_FILE_NAME, "w");
        if (sout == NULL) {
                PRINT_IO_OPEN_ERROR(SUN_POINT_FILE_NAME);
                return;
        }

#ifdef DEBUG
        printf("%s: enrty\n", __func__);
#endif
        unsigned int i;
        double r;
        for (i = 0; i < solution->size; ++i) {
                fprintf(oout, "%lf %lf\n", 
                                storage[i].r, storage[i].theta);
        }

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
        fclose(oout);
        fclose(fout);
}

static double *get_sorted_r(const iteration_storage_t *sorted_st, const size_t size)
{
        double * sorted_r = dv_alloc(sizeof(double) * size);
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
        if (aout == NULL) {
                PRINT_IO_OPEN_ERROR(AVERAGE_R_FILE_NAME);
                return;
        }

        double r = st[0].r;
        int estimate_counter = solution->size;
        int left_bound = 0;
        int i_counter = 1; 

        bool r_into_middle(const double v) {
                return (v < AVERAGE_COUNT_EDGE_R_BOUND) &&
                        (v > AVERAGE_COUNT_EDGE_L_BOUND);
        }

        bool is_last_step(const int c, const int s) {
                return (c - s) <= 0;
        }


        while (estimate_counter > 0) {
                /* Setting for size */
                int size = r_into_middle(r) ? AVERAGE_COUNT_BASE 
                                            : AVERAGE_COUNT_EDGE;

                if (is_last_step(estimate_counter, size)) {
                        size = estimate_counter; 
                }

                /* Get & print */
                a = get_average_theta(&st[left_bound], solution, size);
                double median_r = get_median(&sorted_r[left_bound], size);
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

void dump_background(const iteration_storage_t *st, 
                     const opt_t *solution, 
                     const int b_count)
{       
        FILE *dout = fopen("background.txt", "w");
        if (dout == NULL) {
                PRINT_IO_OPEN_ERROR("background.txt");
                return;
        }

        average_res_t *a;
        double *sorted_r = get_sorted_r(st, solution->size);

        double r = st[0].r;
        int estimate_counter = solution->size;
        int left_bound = 0;
        int i_counter = 1; 

        bool is_last_step(const int c, const int s) {
                return (c - s) <= 0;
        }

        double sd_sd = 0;
        const double sd_total = sqrt(solution->sq / (solution->size + solution->s.size + 1));
        int size = solution->size / b_count; 

        while (estimate_counter > 0) {
                /* Setting for size */

                if (is_last_step(estimate_counter, size)) {
                        size = estimate_counter; 
                }

                /* Get & print */
                a = get_average_theta(&st[left_bound], solution, size);
                sd_sd += pow_double(a->sd - sd_total, 2);

                fprintf(dout, "%lf %lf %lf %lf %d\n",
                               get_median(&sorted_r[left_bound], size),
                               a->err,
                               a->sd,
                               sd_total,
                               size);

                /*  Next step prepare */
                r = sorted_r[left_bound + size - 1];
                left_bound += size;
                estimate_counter -= size;
        }

        FILE *fout = fopen("bk_sd.txt", "w");
        if (dout == NULL) {
                PRINT_IO_OPEN_ERROR("bk_sd.txt");
                return;
        }
        fprintf(fout, "%lf\n", sqrt(sd_sd / (b_count + 1)));

        fclose(fout);
        fclose(dout);
}
