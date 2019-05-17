#include <math.h>
#include <string.h>
#include <stdio.h>
#include "asserts.h"

#include "types.h"
#include "math.h"
#include "opt.h"
#include "unicore.h"
#include "core_b.h"
#include "core_l.h"
#ifdef DEBUG
#include "debug.h"
#endif // DEBUG

void fill_mnk_matrix_vr(linear_equation_t *eq,
                         apogee_rc_table_t *table);
void fill_mnk_matrix(linear_equation_t *eq,
                        apogee_rc_table_t *table, eq_mode_t mode);
/**
 * TODO: generic
 */
static double residuals_line(const linear_eq_solve_t *v,
                             apogee_rc_t *line,
                             const double r_0)
{
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i) {
                mod_v += v->data[i] * core_vr_get_beta_n(line, i);
        }
        for (i = BETA_QTY; i < v->size; ++i) {
                mod_v += core_vr_get_alpha_n(line, i - BETA_QTY + 1, r_0) * v->data[i];
        }

        line->eps = fabs(line->v_helio - mod_v);
        return pow_double(line->v_helio - mod_v, 2);
}

double get_mod_vr(const opt_t *solution,
                  const apogee_rc_t *line)
{
        const double r_0 = GET_SOLUTION_R0(solution);
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i) {
                mod_v += solution->s.data[i] * core_vr_get_beta_n(line, i);
        }
        for (i = BETA_QTY; i < solution->s.size; ++i) {
                mod_v += core_vr_get_alpha_n(line, i - BETA_QTY + 1, r_0) *
                                solution->s.data[i];
        }

        return mod_v;
}

double get_mod_b(const opt_t *solution,
                 const apogee_rc_t *line)
{
        const double r_0 = GET_SOLUTION_R0(solution);

        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i) {
                mod_v += solution->s.data[i] * core_b_get_beta_n(line, i);
        }
        for (i = BETA_QTY; i < solution->s.size; ++i) {
                mod_v += core_b_get_alpha_n(line, i - BETA_QTY + 1, r_0) *
                                solution->s.data[i];
        }

        return mod_v;
}

double get_mod_l(const opt_t *solution,
                 const apogee_rc_t *line)
{
        const double r_0 = GET_SOLUTION_R0(solution);

        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i) {
                mod_v += solution->s.data[i] * core_l_get_beta_n(line, i);
        }
        for (i = BETA_QTY; i < solution->s.size; ++i) {
                mod_v += core_l_get_alpha_n(line, i - BETA_QTY + 1, r_0) *
                                solution->s.data[i];
        }

        return mod_v;
}

static double residuals_summary(const linear_eq_solve_t *v,
                                apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                sum += residuals_line(v, &table->data[i],
                                         GET_TABLE_R0(table));
        }
        assert(sum > 0);
        return sum;
}

double opt_residuals_summary(const linear_eq_solve_t *v,
                             apogee_rc_table_t *table)
{
        return residuals_summary(v, table);
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

#if 0 // TODO: Will be? Now Monte-Karlo give solutions
static inline void step_forward(apogee_rc_table_t *table,
                                const double step, bound_t type)
{
        switch (type) {
                case LOWER:
                        update_table_R0(table, GET_TABLE_R0(table) - step);
                        break;
                case UPPER:
                        update_table_R0(table, GET_TABLE_R0(table) + step);
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
        update_table_R0(table, r_0);

        fill_mnk_matrix_vr(eq, table);
        solve(eq, &s);
        double sq = residuals_summary(&s, table);
        double thr_sq = sq * (1.0 + 1.0 / (table->size + eq->size + 1.0));

        while (step >= SEARCH_PRECISION) {
                step /= STEP_DIVISOR;
                while (1) {
                        step_forward(table, step, type);
                        fill_mnk_matrix_vr(eq, table);
                        solve(eq, &s);
                #ifdef DEBUG
                //        print_vector(s.data, s.size);
                #endif
                        sq = residuals_summary(&s, table);
                #ifdef DEBUG
                        printf("%s: r_0 = %lf, sq = %lf, need >= %lf\n",
                                __func__, GET_TABLE_R0(table), sq, thr_sq);
                #endif
                        if (sq >= thr_sq)
                                break;
                }
                step_backward(table, step, type);
        }

        return GET_TABLE_R0(table);
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
#endif

opt_t *opt_linear(linear_equation_t *eq,
                  apogee_rc_table_t *table,
                  opt_params_t *params)
{
        linear_eq_solve_t s = {
                .data = dv_alloc(sizeof(double) * eq->size),
                .size = eq->size
        };

        assert(params != NULL);
        assert(params->residuals_summary != NULL);
        assert(params->fill_mnk_matrix != NULL);

        double low_r = LOWER_BOUND_R0;
        double high_r = UPPER_BOUND_R0;
        double step = (high_r - low_r) / 16;

        update_table_R0(table, low_r);

        params->fill_mnk_matrix(eq, table);
        solve(eq, &s);
        double sq = params->residuals_summary(&s, table);

        opt_t opt_params = {
                .s = { 0 },
                .r_0 = GET_TABLE_R0(table),
                .sq = sq,
                .size = table->size
        };
        opt_params.s.data = dv_alloc(sizeof(double) * s.size);

        while (step > SEARCH_PRECISION) {
                while (low_r < high_r) {
                        update_table_R0(table, low_r);
                        params->fill_mnk_matrix(eq, table);
                        solve(eq, &s);

                        double sq_tmp = params->residuals_summary(&s, table);
                        if (sq_tmp < opt_params.sq) {
                                opt_params.s = s;
                                opt_params.r_0 = GET_TABLE_R0(table);
                                opt_params.sq = sq_tmp;
                                opt_params.s.data = dv_alloc(sizeof(double) * s.size);
                                memcpy(opt_params.s.data, s.data, s.size * sizeof(double));
                        }

			low_r += step;
                }

                low_r = opt_params.r_0 - step;
                high_r = opt_params.r_0 + step;
                step /= STEP_DIVISOR;
        }


	update_table_R0(table, opt_params.r_0);
        opt_t *ret = dv_alloc(sizeof(opt_t));
        *ret = opt_params;
        return ret;
}


/**
 * @table: data for research
 * @f: optimized function
 * @precalc_errors: criteria for exception
 */
opt_t *exception_algorithm(apogee_rc_table_t *table,
                           opt_t *(*f)(apogee_rc_table_t *),
                           void (*precalc_errors)(apogee_rc_table_t *,
                                                  const double))
{
        parser_t *cfg = get_parser();
        cfg->filter = MATCH_FILTER;
        unsigned int old_size = table->size;
        opt_t *solution = f(table);

        while (true) {
                precalc_errors(table, get_limit_by_eps(table->size));
                filter_get_and_apply(table);
                if (table->size == old_size)
                        break;
                solution = f(table);
                old_size = table->size;
        }

        return solution;
}

/**
 * TODO: need it?
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
*/
