#include <math.h>
#include <assert.h>
#include <string.h>

#include "types.h"
#include "mem.h"
#include "db.h"
#include "core.h"
#include "core_l.h"
#include "graph.h"
#include "generators.h"
#include "trigonometry.h"
#include "unicore.h"
#include "utils.h"

static matrix_line_t g_matrix_line;

double core_l_get_beta_n(const apogee_rc_t *line, beta_ord_t type)
{
        switch (type) {
                case FIRST:
                        return line->sin_l / (line->dist);
                case SECOND:
                        return -(line->cos_l) / (line->dist);
                case THIRD:
                        return -(line->cos_b);
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

        const double tmp = (r0 * line->cos_l / r - (line->cos_b)) / R;

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

                for (i = BETA_QTY; i < len; ++i) {
                        line->_[i] = core_l_get_alpha_n(&table->data[j], i - BETA_QTY + 1, table->r_0);
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
                mod_v += core_l_get_alpha_n(line, i - BETA_QTY + 1, r_0) *
                                solution->s.data[i];
        }

        return mod_v;
}

static double residuals_line(const linear_eq_solve_t *v,
                             apogee_rc_t *line,
                             const double r_0)
{
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i) {
                mod_v += v->data[i] * core_l_get_beta_n(line, i);
        }
        for (i = BETA_QTY; i < v->size; ++i) {
                mod_v += core_l_get_alpha_n(line, i - BETA_QTY + 1, r_0) * v->data[i];
        }

        line->eps = fabs(line->pm_l - mod_v);
        return pow_double(line->pm_l - mod_v, 2);
}

static double residuals_summary_opt(const linear_eq_solve_t *v,
				    apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                sum += residuals_line(v, &table->data[i], table->r_0);
        }
        assert(sum > 0);
        return sum;
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


void precalc_errors_mu_l(apogee_rc_table_t *table,
                         const double limit)
{
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                if (table->data[i].eps / table->sigma[L_PART] > limit) {
                        table->data[i].pm_match = 0;
                }
        }
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

static opt_t *core_l_opt_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        double *matrix = dv_alloc(sizeof(double) * (size + BETA_QTY) *
                                                   (size + BETA_QTY));
        linear_equation_t eq = {
                .data = matrix,
                .right = dv_alloc(sizeof(double) * (size + BETA_QTY)),
                .size = size + BETA_QTY,
                .ord = size
        };

        opt_params_t params = {
                .residuals_summary = residuals_summary_opt,
                .fill_mnk_matrix = core_l_fill_mnk_matrix,
        };
        opt_t *solution = opt_linear(&eq, table, &params);

        table->omega_0 = solution->s.data[BETA_QTY - 1];
        table->r_0 = solution->r_0;
        table->sigma[L_PART] = solution->sq / (table->size - eq.size - 1);
        solution->sq = sqrt(table->sigma[L_PART]);
	return solution;
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
        opt_params.sq = sqrt(opt_params.sq / (table->size - s.size - 1));

        *ret = opt_params;
        return ret;
}

opt_t *core_l_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        cfg->filter = MATCH_FILTER;
        table = get_limited_generic(table, filter_factory(cfg), L_FILTER);

        unsigned int size = cfg->ord;
        double *matrix = dv_alloc(sizeof(double) * (size + BETA_QTY) *
                                                   (size + BETA_QTY));
        linear_equation_t eq = {
                .data = matrix,
                .right = dv_alloc(sizeof(double) * (size + BETA_QTY)),
                .size = size + BETA_QTY,
                .ord = size
        };

        opt_t *opt = core_l_get_linear_solution(&eq, table);
        table->sigma[L_PART] = pow_double(opt->sq, 2);
        table->omega_0 = opt->s.data[BETA_QTY - 1];
#ifdef DEBUG_L
        dump_core_l_solution(opt);
        dump_table_parameters(table, NULL);
#endif
        return opt;
}


void get_partial_l_solution(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        double *matrix = dv_alloc(sizeof(double) * (size + BETA_QTY) *
                                                   (size + BETA_QTY));
        linear_equation_t eq = {
                .data = matrix,
                .right = dv_alloc(sizeof(double) * (size + BETA_QTY)),
                .size = size + BETA_QTY,
                .ord = size
        };

        cfg->filter = MATCH_FILTER;
        table = get_limited_generic(table, filter_factory(cfg), L_FILTER);

        opt_params_t params = {
                .residuals_summary = residuals_summary_opt,
                .fill_mnk_matrix = core_l_fill_mnk_matrix,
        };
        opt_t *solution = opt_linear(&eq, table, &params);

        cfg->filter = ERR_FILTER;
        cfg->h = get_limit_by_eps(table->size);
        cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
        unsigned int old_size = table->size;
        while (true) {
		filter_get_and_apply(table);
                solution = opt_linear(&eq, table, &params);
                cfg->h = get_limit_by_eps(table->size);
                cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
                if (table->size == old_size)
                        break;
                old_size = table->size;
        }

        table->omega_0 = solution->s.data[BETA_QTY - 1];

        table->r_0 = solution->r_0;
        table->sigma[L_PART] = solution->sq / (table->size - eq.size - 1);
        solution->sq = sqrt(table->sigma[L_PART]);

        dump_rotation_curve_l(solution);
        dump_objects_theta_R(table, solution, L_PART, "l_objs.txt");

        if (cfg->draw_profile)
                dump_profile(&eq, table, &params, "l_profile.txt");

        mk_params_t mk_params = {
                .f_entry = core_l_opt_entry,
                .f_point_by_solution = get_point_by_l_solution,
                .f_table_by_solution = fill_table_by_l_solution,
                .count = cfg->mksize,
        };

        opt_t *mk_sol = monte_carlo_entry(solution,
                                          table,
                                          &mk_params);


        dump_result(mk_sol);
        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_theta_R(dumped, solution, L_PART, "l_objs_err.txt");
        dump_vr_solution(mk_sol);
        dump_objects_xyz(dumped, dumped->size, "ERROR_LIMITED");
}

