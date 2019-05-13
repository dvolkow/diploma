#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#include "db.h"
#include "core.h"
#include "core_b.h"
#include "core_l.h"
#include "core_vr.h"
#include "types.h"
#include "io.h"
#include "mem.h"
#include "math.h"
#include "debug.h"
#include "opt.h"
#include "graph.h"
#include "utils.h"
#include "unicore.h"
#include "generators.h"

static matrix_line_t g_matrix_line;

double core_vr_get_beta_n(const apogee_rc_t *line, beta_ord_t type)
{
        switch (type) {
                case FIRST:
                        return -line->cos_l * line->cos_b;
                case SECOND:
                        return -line->sin_l * line->cos_b;
                case THIRD:
                        return -line->sin_b;
                default:
#ifdef DEBUG
                        printf("%s: type error!\n", __func__);
#endif
                        return 0;
        }
}

static inline double s_alpha_n(const double R,
                               const double sinl,
                               const double cosb,
                               const double r_0, const unsigned int n)
{
        if (n == 1)
                return -2 * (R - r_0) * r_0 * sinl * cosb / R;
        else 
                return r_0 * pow_double(R - r_0, n) * sinl * cosb / (R * dv_factorial(n));
}

double core_vr_get_alpha_n(const apogee_rc_t *line,
                           const unsigned int n,
                           const double r_0)
{
        assert(n > 0);

        double R = get_R_distance(line, r_0);
        return s_alpha_n(R, line->sin_l, line->cos_b, r_0, n);
}

void fill_mnk_matrix_vr(linear_equation_t *eq,
                        apogee_rc_table_t *table)
{
        unsigned int i, j, k;
        unsigned int len = eq->size;
        matrix_line_t *line = &g_matrix_line;

        assert(len <= BETA_QTY + MAX_ORDER_SOLUTION);

        memset(eq->data, 0, sizeof(double) * eq->size * eq->size);
        memset(eq->right, 0, sizeof(double) * eq->size);

        for (j = 0; j < table->size; ++j) {
                for (i = 0; i < BETA_QTY; ++i) {
                        line->_[i] = core_vr_get_beta_n(&table->data[j], i);
                }

                for (i = BETA_QTY; i < eq->size; ++i) {
                        line->_[i] = core_vr_get_alpha_n(&table->data[j], i - BETA_QTY + 1, table->r_0);
                }

                for (i = 0; i < len; ++i) {
                        double m = line->_[i];
                        for (k = 0; k < len; ++k) {
                                eq->data[i * len + k] += line->_[k] * m;
                        }
                        eq->right[i] += table->data[j].v_helio * m;
                }
        }
#ifdef DEBUG
        print_matrix(eq->data, len);
#endif
}

void fill_mnk_matrix(linear_equation_t *eq,
                        apogee_rc_table_t *table, eq_mode_t mode)
{
        switch (mode) {
                case VR_MODE:
                        fill_mnk_matrix_vr(eq, table);
                        break;
                default:
                        printf("%s: mode not implemented!\n", __func__);
                        exit(1);
        }
}

void filter_get_and_apply(apogee_rc_table_t *table)
{
        /* TODO: Filter assigner: */
        parser_t *cfg = get_parser();

        switch (cfg->filter) {
                case L_FILTER:
                case B_FILTER:
                case ERR_FILTER:
                case MATCH_FILTER:
                        table = get_limited_generic(table,
                                                    filter_factory(cfg),
                                                    L_FILTER);
                        break;
                default:
                        return;
        }
}

static opt_t *__get_solution_iterate(apogee_rc_table_t *table,
                                        linear_equation_t *eq)
{
        assert(table->size != 0);

        filter_get_and_apply(table);

        opt_params_t params = {
                .residuals_summary = opt_residuals_summary,
                .fill_mnk_matrix = fill_mnk_matrix_vr,
        };

        opt_t *solution = opt_linear(eq, table, &params);
        return solution;
}

void get_iterate_solution(apogee_rc_table_t *table,
                          opt_t *solution);
void get_iterate_solution_nerr(apogee_rc_table_t *table,
                               opt_t *solution);

static opt_t *vr_partial_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        unsigned int dim = size + BETA_QTY;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_params_t params = {
                .residuals_summary = opt_residuals_summary,
                .fill_mnk_matrix = fill_mnk_matrix_vr,
        };

        opt_t *opt = opt_linear(&eq, table, &params);
        table->r_0 = opt->r_0;
        opt->sq = sqrt(opt->sq / (table->size - opt->s.size - 1));

        return opt;

}

void get_partial_vr_solution(apogee_rc_table_t *table)
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
                .residuals_summary = opt_residuals_summary,
                .fill_mnk_matrix = fill_mnk_matrix_vr,
        };
        opt_t *solution = opt_linear(&eq, table, &params);

        cfg->filter = ERR_FILTER;
        cfg->h = get_limit_by_eps(table->size);
        cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
        unsigned int old_size = table->size;
        while (true) {
                solution = __get_solution_iterate(table, &eq);
                cfg->h = get_limit_by_eps(table->size);
                cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
                if (table->size == old_size)
                        break;
                old_size = table->size;
        }

        table->omega_0 = OMEGA_SUN - solution->s.data[V_P] / solution->r_0;

        table->r_0 = solution->r_0;
        table->sigma[VR_PART] = solution->sq / (table->size - eq.size - 1);
        solution->sq = sqrt(table->sigma[VR_PART]);
        printf("%s: sq = %lf\n", __func__, solution->sq);

        dump_rotation_curve_vr(solution);
        dump_objects_theta_R(table, solution, TOTAL_QTY, "vr_objs.txt");

        mk_params_t mk_params = {
                .f_entry = vr_partial_entry,
                .f_point_by_solution = get_point_by_vr_solution,
                .f_table_by_solution = fill_table_by_vr_solution,
                .count = cfg->mksize,
        };

        opt_t *mk_sol = monte_carlo_entry(solution,
                                          table,
                                          &mk_params);


        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_theta_R(dumped, solution, TOTAL_QTY, "vr_objs_err.txt");
        dump_vr_solution(mk_sol);
        dump_objects_xyz(dumped, dumped->size, "ERROR_LIMITED");
}

void get_united_solution(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        double *matrix = dv_alloc(sizeof(double) * (size + BETA_QTY) *
                                                   (size + BETA_QTY));

        filter_get_and_apply(table);

        dump_objects_xyz(table, table->size, name_for_obj(178, 0, __func__));
        dump_table(table);

        linear_equation_t eq = {
                .data = matrix,
                .right = dv_alloc(sizeof(double) * (size + BETA_QTY)),
                .size = size + BETA_QTY,
                .ord = size
        };

        opt_params_t params = {
                .residuals_summary = opt_residuals_summary,
                .fill_mnk_matrix = fill_mnk_matrix_vr,
        };
        opt_t *solution = opt_linear(&eq, table, &params);
        // TODO: function's interface need improve
        cfg->filter = ERR_FILTER;
        cfg->h = get_limit_by_eps(table->size);
        cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
        unsigned int old_size = table->size;
        while (true) {
                solution = __get_solution_iterate(table, &eq);
                cfg->h = get_limit_by_eps(table->size);
                cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
                if (table->size == old_size)
                        break;
                old_size = table->size;
        }

        get_iterate_solution(table, solution);
        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_xyz(dumped, dumped->size, "ERROR_LIMITED");
        //get_iterate_solution_nerr(table, solution);
}


#if 0
void get_iterate_solution_nerr(apogee_rc_table_t *table,
                               opt_t *solution)
{
        parser_t *cfg = get_parser();
        cfg->filter = MATCH_FILTER;
        table = get_limited_generic(table, filter_factory(cfg), L_FILTER);

        solution = united_with_nature_errs_entry(table);

        cfg->filter = ERR_FILTER;
        cfg->h = get_limit_by_eps(table->size);
        cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
        unsigned int old_size = table->size;
        unsigned int i = 1;
        while (true) {
                printf("%s: iteration #%u, size %u\n", __func__, i++, table->size);
                filter_get_and_apply(table);
                if (table->size == old_size)
                        break;
                solution = united_with_nature_errs_entry(table);
                cfg->h = get_limit_by_eps(table->size);
                cfg->l = sqrt(solution->sq / (solution->size - solution->s.size - 1));
                old_size = table->size;
        }
        
        mk_params_t mk_params = {
                .f_entry = united_with_nature_errs_entry,
                .f_point_by_solution = get_point_by_uni_solution,
                .f_table_by_solution = fill_table_by_uni_solution,
                .count = cfg->mksize,
        };

        opt_t *mk_sol = monte_carlo_entry(solution,
                                          table,
                                          &mk_params);
        dump_united_solution(mk_sol);
}
#endif


void vr_b_iterations(apogee_rc_table_t *table)
{
        double sd[TOTAL_QTY];
        table->w_sun = W_SUN_START;
        double w_old = W_SUN_START;
        double r_old = table->r_0;

        core_vr_entry(table);
        double r_new = table->r_0;
        core_b_entry(table);
        double w_new = table->w_sun;
        sd[B_PART] = table->sigma[B_PART];

        while (ITER_CONDITION(w_old, w_new, r_old, r_new)) {
                r_old = r_new;
                w_old = w_new;
                core_vr_entry(table);
                r_new = table->r_0;
                core_b_entry(table);
                w_new = table->w_sun;
                sd[B_PART] = table->sigma[B_PART];
        }

        core_vr_entry(table);
        sd[VR_PART] = table->sigma[VR_PART];
        core_l_entry(table);
        sd[L_PART] = table->sigma[L_PART];
        uni_g_sd_init(sd);
}


void get_iterate_solution(apogee_rc_table_t *table,
                          opt_t *solution)
{
        double sd[TOTAL_QTY];

        parser_t *cfg = get_parser();
        cfg->filter = MATCH_FILTER;
        table = get_limited_generic(table, filter_factory(cfg), L_FILTER);

        table->w_sun = W_SUN_START;
        double w_old = W_SUN_START;
        double r_old = table->r_0;
        solution = exception_algorithm(table,
                                       core_vr_entry,
                                       precalc_errors_vr);
        sd[VR_PART] = table->sigma[VR_PART];
        double r_new = solution->r_0;
        solution = exception_algorithm(table,
                                       core_b_entry,
                                       precalc_errors_mu_b);
        double w_new = table->w_sun;
        sd[B_PART] = table->sigma[B_PART];

        while (ITER_CONDITION(w_old, w_new, r_old, r_new)) {
                r_old = r_new;
                w_old = w_new;
                solution = exception_algorithm(table,
                                               core_vr_entry,
                                               precalc_errors_vr);
                r_new = table->r_0;
                sd[VR_PART] = table->sigma[VR_PART];
                solution = exception_algorithm(table,
                                               core_b_entry,
                                               precalc_errors_mu_b);
                sd[B_PART] = table->sigma[B_PART];
                w_new = table->w_sun;
        }

        solution = core_vr_entry(table);
        sd[VR_PART] = table->sigma[VR_PART];

        solution = core_b_entry(table);
        sd[B_PART] = table->sigma[B_PART];

        solution = core_l_entry(table);
        sd[L_PART] = table->sigma[L_PART];
        dump_objects_theta_R(table, solution, L_PART, "L_PART_OBJ.txt");
        dump_part_rotation_curve(solution, L_PART, "l_cur.txt", table->omega_0);


        solution = core_vr_entry(table);
        dump_objects_theta_R(table, solution, VR_PART, "VR_PART_OBJ.txt");
        dump_part_rotation_curve(solution, VR_PART, "vr_cur.txt", table->omega_0);

        solution = core_b_entry(table);
        dump_objects_theta_R(table, solution, B_PART, "B_PART_OBJ.txt");
        dump_part_rotation_curve(solution, B_PART, "b_cur.txt", table->omega_0);

        uni_g_sd_init(sd);
        // TODO: need modify dispersion?
        if (cfg->bolter > 0)
                solution = exception_algorithm(table,
                                               united_entry,
                                               precalc_errors_uni); 
        else
                solution = united_entry(table);
 
        mk_params_t mk_params = {
                .f_entry = united_entry,
                .f_point_by_solution = get_point_by_uni_solution,
                .f_table_by_solution = fill_table_by_uni_solution,
                .count = cfg->mksize,
        };

        opt_t *mk_sol = monte_carlo_entry(solution,
                                          table,
                                          &mk_params);
        dump_united_solution(mk_sol);
        dump_table_parameters(table, mk_sol);

        dump_result(mk_sol);
        dump_uni_rotation_objs(table, solution);

        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_uni_rotation_objs_named(dumped,
                                     solution,
                                     "ERROR_LIMITED_R_THETA");
}
