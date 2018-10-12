#include <math.h>
#include <string.h>
#include "asserts.h"

#include "types.h"
#include "math.h"
#include "opt.h"
#ifdef DEBUG
#include "debug.h"
#endif // DEBUG


double get_beta_n(const apogee_rc_t *line, beta_ord_t type);
double get_alpha_n(const apogee_rc_t *line, 
                   const double r_0,
                   const int n);

void fill_mnk_matrix_vr(linear_equation_t *eq,
                         apogee_rc_table_t *table);
void fill_mnk_matrix(linear_equation_t *eq, 
                        apogee_rc_table_t *table, eq_mode_t mode);
/**
 * TODO: generic
 */
static double residuals_line(const linear_equation_t *eq,
                             const linear_eq_solve_t *v,
                             apogee_rc_t *line,
                             const double r_0)
{
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i) {
                mod_v += v->data[i] * get_beta_n(line, i);
        }
        for (i = BETA_QTY; i < v->size; ++i) {
                mod_v += get_alpha_n(line, r_0, i - BETA_QTY + 1) * v->data[i];
        }

        line->eps = fabs(line->v_helio - mod_v);
        return pow_double(line->v_helio - mod_v, 2);
}

static double residuals_summary(const linear_equation_t *eq, 
                                const linear_eq_solve_t *v, 
                                const apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                sum += residuals_line(eq, v, &table->data[i], table->r_0);
        }
        assert(sum > 0);
        return sum;
}

static void fill_estimated_vector(const apogee_rc_table_t *table,
                                linear_eq_solve_t *v)
{
        unsigned int i;
        assert(table->size == v->size);
        for (i = 0; i < table->size; ++i) {
                v->data[i] = table->data[i].v_helio;
        }
}

static inline void step_forward(apogee_rc_table_t *table,
                                const double step, bound_t type)
{
        switch (type) {
                case LOWER:
                        table->r_0 -= step;
                        break;
                case UPPER:
                        table->r_0 += step;
                        break;
                default:
                        printf("%s: error bound type!\n",
                                        __func__);
        }
}

static inline void step_backward(apogee_rc_table_t *table,
                                const double step, bound_t type)
{
        step_forward(table, step, type == LOWER ? UPPER 
                                                : LOWER);
}

static double __bound_parameter(linear_equation_t *eq, 
                                        apogee_rc_table_t *table,
                                        double r_0, bound_t type) 
{
        linear_eq_solve_t s = {
                .data = dv_alloc(sizeof(double) * eq->size),
                .size = eq->size
        };

        double step = (UPPER_BOUND_R0 - LOWER_BOUND_R0) / 16;
        table->r_0 = r_0;

        fill_mnk_matrix_vr(eq, table);
        solve(eq, &s); 
        double sq = residuals_summary(eq, &s, table);
        double thr_sq = sq * (1.0 + 1.0 / (table->size + eq->size + 1.0));

        while (step >= SEARCH_PRECISION) {
                step /= STEP_DIVISOR;
                while (1) {
                        step_forward(table, step, type);
                        fill_mnk_matrix_vr(eq, table);
                        solve(eq, &s);
                #ifdef DEBUG
                        print_vector(s.data, s.size);
                #endif
                        sq = residuals_summary(eq, &s, table);
                #ifdef DEBUG
                        printf("%s: r_0 = %lf, sq = %lf, need >= %lf\n",
                                __func__, table->r_0, sq, thr_sq);
                #endif
                        if (sq >= thr_sq)
                                break;
                }
                step_backward(table, step, type);
        }

        return table->r_0;
}

double lower_bound_search(linear_equation_t *eq, 
                                apogee_rc_table_t *table,
                                double r_0)
{
        return __bound_parameter(eq, table, r_0, LOWER);
}

double upper_bound_search(linear_equation_t *eq, 
                                apogee_rc_table_t *table,
                                double r_0)
{
        return __bound_parameter(eq, table, r_0, UPPER);
}

opt_t *opt_linear(linear_equation_t *eq, 
                                apogee_rc_table_t *table) 
{
        linear_eq_solve_t s = {
                .data = dv_alloc(sizeof(double) * eq->size),
                .size = eq->size
        };

        double low_r = LOWER_BOUND_R0;
        double high_r = UPPER_BOUND_R0;
        double step = (high_r - low_r) / 32;

        table->r_0 = low_r;

        fill_mnk_matrix_vr(eq, table);
        solve(eq, &s); 
#ifdef DEBUG
        print_vector(s.data, s.size);
#endif
        double sq = residuals_summary(eq, &s, table);

        opt_t opt_params = {
                .s = dv_alloc(sizeof(linear_eq_solve_t)),
                .r_0 = table->r_0,
                .sq = sq,
                .size = table->size
        };
        opt_params.s.data = dv_alloc(sizeof(double) * s.size);

        while (step > SEARCH_PRECISION) {
                while (table->r_0 < high_r) {
                        table->r_0 += step;
                        fill_mnk_matrix_vr(eq, table);
                        solve(eq, &s);

        #ifdef DEBUG
                        print_vector(s.data, s.size);
        #endif
                        double sq_tmp = residuals_summary(eq, &s, table);
        #ifdef DEBUG
                        printf("%s: sd = %lf, r_0 = %lf\n", 
                                        __func__, sq_tmp, table->r_0);
        #endif
                        if (sq_tmp < opt_params.sq) {
                                opt_params.s = s,
                                opt_params.r_0 = table->r_0;
                                opt_params.sq = sq_tmp;
                                opt_params.s.data = dv_alloc(sizeof(double) * s.size);
                                memcpy(opt_params.s.data, s.data, s.size * sizeof(double));
                        }
                }

                low_r = opt_params.r_0 - step;
                high_r = opt_params.r_0 + step;
                step /= STEP_DIVISOR;
                table->r_0 = low_r;
        }


        opt_t *ret = dv_alloc(sizeof(opt_t));
        *ret = opt_params;
        return ret;
}

void get_errors(opt_t *solution, apogee_rc_table_t *table) 
{
        linear_equation_t m, invm;
        m.data = dv_alloc(sizeof(double) * solution->s.size * solution->s.size);
        m.right = dv_alloc(sizeof(double) * solution->s.size);
        m.size = solution->s.size;
        invm.data = dv_alloc(sizeof(double) * solution->s.size * solution->s.size);
        invm.size = solution->s.size;

        fill_mnk_matrix(&m, table, VR_MODE);
        inverse_and_diag(&m, &invm);

        solution->bounds = dv_alloc(sizeof(prec_t) * solution->s.size);
        unsigned int i;
        for (i = 0; i < solution->s.size; ++i) {
                solution->bounds[i].l = 
                        get_error_mnk_estimated(invm.data[i * invm.size + i], invm.size + table->size + 1, solution->sq / (table->size + 1));
        }
}
