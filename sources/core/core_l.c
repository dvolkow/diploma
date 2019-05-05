#include <math.h>
#include <assert.h>
#include "types.h"
#include "mem.h"
#include "core_l.h"
#include "graph.h"

static matrix_line_t g_matrix_line;

double core_l_get_beta_n(const apogee_rc_t *line, beta_ord_t type)
{
        switch (type) {
                case FIRST:
                        return sin(line->l) / (line->dist);
                case SECOND:
                        return -cos(line->l) / (line->dist);
                case THIRD:
                        return -cos(line->b);
                default:
#ifdef DEBUG
                        printf("%s: type error!\n", __func__);
#endif
                        return 0;
        }
}

double core_l_get_alpha_n(const apogee_rc_t *line,
                          const unsigned int n,
                          const double r0)
{
        const double R = get_R_distance(line, r0);
        const double r = line->dist;
        const double delta_R = R - r0;

        const double tmp = (r0 * cos(line->l) / r - cos(line->b)) / R;

        if (n > 1) {
                return (pow_double(delta_R, n) * tmp) / dv_factorial(n);
        } else {
                return -2 * delta_R * tmp;
        }
}

void core_l_fill_mnk_matrix(linear_equation_t *eq,
                            apogee_rc_table_t *table)
{
        unsigned int i, j, k;
        unsigned int len = eq->size;
        matrix_line_t *line = &g_matrix_line;

        memset(eq->data, 0, sizeof(double) * eq->size * eq->size);
        memset(eq->right, 0, sizeof(double) * eq->size);

        for (j = 0; j < table->size; ++j) {
                for (i = 0; i < BETA_QTY; ++i) {
                        line->_[i] = core_l_get_beta_n(&table->data[j], i);
                }

                for (i = BETA_QTY; i < eq->size; ++i) {
                        line->_[i] = core_l_get_alpha_n(&table->data[j], table->r_0, i - BETA_QTY + 1);
                }

                for (i = 0; i < len; ++i) {
                        double m = line->_[i];
                        for (k = 0; k < len; ++k) {
                                eq->data[i * len + k] += line->_[k] * m;
                        }
                        eq->right[i] += table->data[j].pm_l * m;
                }
        }
}


static double core_l_get_mod_v(const opt_t *solution,
                               const apogee_rc_t *line)
{
        const double r_0 = solution->r_0;
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i) {
                mod_v += solution->s.data[i] * core_l_get_beta_n(line, i);
        }

        for (i = BETA_QTY; i < solution->s.size; ++i) {
                mod_v += core_l_get_alpha_n(line, r_0, i - BETA_QTY + 1) * 
                                solution->s.data[i];
        }

        return mod_v;
}


static double core_l_residuals_line(const opt_t *solution,
                                    apogee_rc_t *line)
{
        const double mod_v = core_l_get_mod_v(solution, line);
        /* pm_l and pm_b already multiplied in k */
        line->eps = fabs(line->pm_l - mod_v);
        return pow_double(line->pm_l - mod_v, 2);
}


static double residuals_summary(const opt_t *solution, 
                                const apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                sum += core_l_residuals_line(solution, &table->data[i]);
        }
        assert(sum > 0);
        return sum;
}

void core_l_get_errors(opt_t *solution, apogee_rc_table_t *table)
{
        linear_equation_t m, invm;
        m.data = dv_alloc(sizeof(double) * solution->s.size * solution->s.size);
        m.right = dv_alloc(sizeof(double) * solution->s.size);
        m.size = solution->s.size;
        invm.data = dv_alloc(sizeof(double) * solution->s.size * solution->s.size);
        invm.size = solution->s.size;

        core_l_fill_mnk_matrix(&m, table);
        inverse_and_diag(&m, &invm);

        solution->bounds = dv_alloc(sizeof(prec_t) * solution->s.size);
        unsigned int i;
        for (i = 0; i < solution->s.size; ++i) {
                solution->bounds[i].l =
                        get_error_mnk_estimated(invm.data[i * invm.size + i], invm.size + table->size + 1, solution->sq / (table->size + 1));
        }
}


opt_t *core_l_get_linear_solution(linear_equation_t *eq,
                                  apogee_rc_table_t *table)
{
        opt_t *ret = dv_alloc(sizeof(opt_t));
        linear_eq_solve_t s = {
                .data = dv_alloc(sizeof(double) * eq->size),
                .size = eq->size
        };

        core_l_fill_mnk_matrix(eq, table);
        solve(eq, &s);

        opt_t opt_params = {
                .s = s,
                .r_0 = table->r_0,
                .size = table->size
        };
        opt_params.sq = residuals_summary(&opt_params, table);
        core_l_get_errors(&opt_params, table);
        /* To dump get sd: */
        opt_params.sq = sqrt(opt_params.sq / (table->size + s.size - 1));

        *ret = opt_params;
        return ret;
}

opt_t *core_l_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
//       cfg->filter = MATCH_FILTER;
//       table = get_limited_generic(table, filter_factory(cfg), L_FILTER);

        int size = cfg->ord;
        double *matrix = dv_alloc(sizeof(double) * (size + BETA_QTY) *
                                                   (size + BETA_QTY));
        linear_equation_t eq = {
                .data = matrix,
                .right = dv_alloc(sizeof(double) * (size + BETA_QTY)),
                .size = size + BETA_QTY,
                .ord = size
        };

        opt_t *opt = core_l_get_linear_solution(&eq, table);
        dump_core_l_solution(opt);
        return opt;
}