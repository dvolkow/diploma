#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#include "db.h"
#include "debug.h"
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
                        return 0;
        }
}

static inline double s_alpha_n(const double R,
                               const double sinl,
                               const double cosb,
                               const double r_0, const unsigned int n)
{
        if (n > 1)
                return r_0 * pow_double(R - r_0, n) * sinl * cosb / (R * dv_factorial(n));
        else
                return -2 * (R - r_0) * r_0 * sinl * cosb / R;
}

double core_vr_get_alpha_n(const apogee_rc_t *line,
                           const unsigned int n,
                           const double r_0)
{
        assert(n > 0);

#ifndef PRECACHED_TABLE_R
        const double R = get_R_distance(line, r_0);
#else
        const double R = GET_LINE_R(line);
#endif
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
                        line->_[i] = core_vr_get_alpha_n(&table->data[j], i - BETA_QTY + 1, GET_TABLE_R0(table));
                }

                for (i = 0; i < len; ++i) {
                        double m = line->_[i];
                        for (k = 0; k < len; ++k) {
                                eq->data[i * len + k] += line->_[k] * m;
                        }
                        eq->right[i] += table->data[j].v_helio * m;
                }
        }
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

void get_iterate_solution(apogee_rc_table_t *table);
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
        update_table_R0(table, GET_SOLUTION_R0(opt));
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

        table->omega_0 = OMEGA_SUN - solution->s.data[V_P] / GET_SOLUTION_R0(solution);

        update_table_R0(table, GET_SOLUTION_R0(solution));
        table->sigma[VR_PART] = solution->sq / (table->size - eq.size - 1);
        solution->sq = sqrt(table->sigma[VR_PART]);
        printf("%s: sq = %lf\n", __func__, solution->sq);

        dump_rotation_curve_vr(solution);
        dump_objects_theta_R(table, solution, TOTAL_QTY, "vr_objs.txt");

        if (cfg->draw_profile)
                dump_profile(&eq, table, &params, "vr_profile.txt");

        mk_params_t mk_params = {
                .f_entry = vr_partial_entry,
                .f_point_by_solution = get_point_by_vr_solution,
                .f_table_by_solution = fill_table_by_vr_solution,
                .count = cfg->mksize,
		.mul_unfres_name = "mk_results.txt"
        };

        opt_t *mk_sol = monte_carlo_entry(solution,
                                          table,
                                          &mk_params);


        if (mk_sol) dump_result(mk_sol);

        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_theta_R(dumped, solution, TOTAL_QTY, "vr_objs_err.txt");
        if (mk_sol) dump_vr_solution(mk_sol);
        dump_objects_xyz(dumped, dumped->size, "ERROR_LIMITED");

	find_r_0_bounds(table, solution, &params, &eq);
	dump_r0_bounds(solution);
}


void get_united_solution(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();

        cfg->filter = MATCH_FILTER;
        table = get_limited_generic(table, filter_factory(cfg), L_FILTER);
        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_xyz(dumped, dumped->size, "missing_xyz.txt");
	db_clear(ERROR_LIMITED);

        get_iterate_solution(table);
}

void find_united_sigma_0_solution(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        cfg->filter = MATCH_FILTER;
        table = get_limited_generic(table, filter_factory(cfg), L_FILTER);

        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_xyz(dumped, dumped->size, "missing_xyz.txt");
	db_clear(ERROR_LIMITED);


        if (cfg->bolter > 0)
                vr_b_iterations(table);

	double sigma_0_opt = cfg->sigma_0;
        opt_t *solution = united_with_nature_errs_entry(table);

	double chi_opt = solution->sq;
	// warning: not -b option here!
	const double n_free = 3 * table->size - solution->s.size - 1;
	double step = 1;

	double start = cfg->sigma_0_l;
	double end = cfg->sigma_0_h;

	opt_t opt_solution = *solution;
	double diff_opt = fabs(chi_opt - n_free);

	while (step > SEARCH_PRECISION) {
		while(start < end) {
			cfg->sigma_0 = start;
			solution = united_with_nature_errs_entry(table);
			double diff = fabs(solution->sq - n_free);
			if (diff < diff_opt) {
				printf("%s: old diff %lf, new diff %lf\n",
				       __func__, diff_opt, diff);
				diff_opt = diff;
				chi_opt = solution->sq;
				opt_solution = *solution;
				sigma_0_opt = cfg->sigma_0;
			}
			start += step;
		}

		start = sigma_0_opt - step;
		end = sigma_0_opt + step;
		step /= STEP_DIVISOR;
	}

	printf("%s: optimal sigma_0 %lf\n", __func__, sigma_0_opt);
	partial_dump_unfriendly_result(&opt_solution, MAX_PRINTF_COLS + 1, "u_unfresult.txt");
	precalc_vsd_to_dump(table);
	dump_residuals(table);
        dump_uni_rotation_objs(table, solution);

        mk_params_t mk_params = {
                .f_entry = united_with_nature_errs_entry,
                .f_point_by_solution = get_point_by_uni_solution,
                .f_table_by_solution = fill_table_by_uni_solution_sigma_0,
                .count = cfg->mksize,
		.mul_unfres_name = "mk_results.txt"
        };

        opt_t *mk_sol = monte_carlo_entry(&opt_solution,
                                          table,
                                          &mk_params);

        if (mk_sol) {
		dump_mk_errors_uni(mk_sol);
		dump_mk_values(mk_sol);
		dump_result(mk_sol);
		dump_united_solution(mk_sol);
	}

        if (cfg->draw_profile) {
                dump_united_sigma0_solution_profile(table, cfg->ord);
        }

	dump_united_solution_r0_bounds(table, solution);
}


void get_united_sigma_0_solution(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        cfg->filter = MATCH_FILTER;
        table = get_limited_generic(table, filter_factory(cfg), L_FILTER);

        opt_t *solution = united_with_nature_errs_entry(table);

	if (cfg->bolter > 0)
                solution = exception_algorithm(table,
                                               united_with_nature_errs_entry,
					       precalc_errors_uni_sigma_0);

	partial_dump_unfriendly_result(solution, MAX_PRINTF_COLS + 1, "u_unfresult.txt");
        dump_uni_rotation_objs(table, solution);

	precalc_vsd_to_dump(table);
	dump_residuals(table);
        if (cfg->draw_profile) {
                dump_united_sigma0_solution_profile(table, cfg->ord);
        }

        mk_params_t mk_params = {
                .f_entry = united_with_nature_errs_entry,
                .f_point_by_solution = get_point_by_uni_solution,
                .f_table_by_solution = fill_table_by_uni_solution_sigma_0,
                .count = cfg->mksize,
		.mul_unfres_name = "mk_results.txt"
        };

        opt_t *mk_sol = monte_carlo_entry(solution,
                                          table,
                                          &mk_params);
	if (mk_sol) {
		dump_result(mk_sol);
		dump_united_solution(mk_sol);
		dump_mk_errors_uni(mk_sol);
		dump_mk_values(mk_sol);
	}
        apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_xyz(dumped, dumped->size, "missing_xyz.txt");
	dump_united_solution_r0_bounds(table, solution);
}


void vr_b_iterations(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        double sd[TOTAL_QTY];

        table->w_sun = W_SUN_START;
        double w_old = W_SUN_START;
        double r_old = GET_TABLE_R0(table);

        opt_t *solution = core_vr_entry(table);
        sd[VR_PART] = solution->sq;
        double r_new = GET_TABLE_R0(table);
        solution = core_b_entry(table);
        sd[B_PART] = solution->sq;
        double w_new = table->w_sun;
        solution = core_l_entry(table);
        sd[L_PART] = solution->sq;
        uni_g_sd_init(sd);


        unsigned int new_table_size = table->size;
        unsigned int old_table_size;
        unsigned int i = 0;
#define MAX_ITERATIONS_COUNT 5
        while (true) {
                old_table_size = new_table_size;
                while (ITER_CONDITION(w_old, w_new, r_old, r_new) && (i < MAX_ITERATIONS_COUNT)) {
                        r_old = r_new;
                        w_old = w_new;
                        solution = cfg->bolter == 0 ? core_vr_entry(table)
                                                    : exception_algorithm(table,
                                                                          core_vr_entry,
                                                                          precalc_errors_vr);
                        r_new = GET_TABLE_R0(table);
                        solution = cfg->bolter == 0 ? core_b_entry(table)
                                                    : exception_algorithm(table,
                                                                          core_b_entry,
                                                                          precalc_errors_mu_b);
                        w_new = table->w_sun;
                        ++i;
                }


                solution = cfg->bolter == 0 ? core_l_entry(table)
                                            : exception_algorithm(table,
                                                                  core_l_entry,
                                                                  precalc_errors_mu_l);

                new_table_size = table->size;

                if (new_table_size == old_table_size) {
                        break;
                }
        }

        sd[B_PART] = table->sigma[B_PART];
        dump_objects_theta_R(table, solution, B_PART, "B_PART_OBJ.txt");
        dump_part_rotation_curve(solution, B_PART, "b_cur.txt", table->omega_0);
	partial_dump_unfriendly_result(solution, MAX_PRINTF_COLS, "b_unfresult.txt");

        solution = core_vr_entry(table);
        sd[VR_PART] = table->sigma[VR_PART];
        dump_objects_theta_R(table, solution, VR_PART, "VR_PART_OBJ.txt");
        dump_part_rotation_curve(solution, VR_PART, "vr_cur.txt", table->omega_0);
	partial_dump_unfriendly_result(solution, MAX_PRINTF_COLS - 1, "vr_unfresult.txt");

        solution = core_l_entry(table);
        sd[L_PART] = table->sigma[L_PART];
        dump_objects_theta_R(table, solution, L_PART, "L_PART_OBJ.txt");
        dump_part_rotation_curve(solution, L_PART, "l_cur.txt", table->omega_0);
	partial_dump_unfriendly_result(solution, MAX_PRINTF_COLS, "l_unfresult.txt");
        uni_g_sd_init(sd);
}


void get_iterate_solution(apogee_rc_table_t *table)
{
	opt_t *solution;
        parser_t *cfg = get_parser();


        vr_b_iterations(table);
        // TODO: need modify dispersion?
        if (cfg->bolter > 0)
                solution = exception_algorithm(table,
                                               united_entry,
                                               precalc_errors_uni);
        else
                solution = united_entry(table);

	precalc_vsd_to_dump(table);
	dump_residuals(table);
        dump_uni_rotation_objs(table, solution); // must be here wher R is precached!

	partial_dump_unfriendly_result(solution, MAX_PRINTF_COLS + 1, "u_unfresult.txt");
        if (cfg->draw_profile)
                dump_united_solution_profile(table, cfg->ord);

        mk_params_t mk_params = {
                .f_entry = united_entry,
                .f_point_by_solution = get_point_by_uni_solution,
                .f_table_by_solution = fill_table_by_uni_solution,
                .count = cfg->mksize,
		.mul_unfres_name = "mk_results.txt"
        };

        opt_t *mk_sol = monte_carlo_entry(solution,
                                          table,
                                          &mk_params);
	if (mk_sol) {
		dump_united_solution(mk_sol);
		dump_table_parameters(table, mk_sol);
		dump_result(mk_sol);
		dump_mk_errors_uni(mk_sol);
		dump_mk_values(mk_sol);
	}

	dump_united_solution_r0_bounds(table, solution);
}
