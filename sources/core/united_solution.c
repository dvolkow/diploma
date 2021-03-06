#include <math.h>
#include <assert.h>

#include "debug.h"
#include "types.h"
#include "core.h"
#include "mem.h"
#include "core_l.h"
#include "core_b.h"
#include "unicore.h"
#include "trigonometry.h"

static matrix_line_t g_matrix_line;
static double g_sq[TOTAL_QTY];
//static linear_equation_t g_tmp_eq;


static double united_l_get_beta_n(const apogee_rc_t *line, beta_ord_t type)
{
        switch (type) {
                case FIRST:
                case SECOND:
                        return core_l_get_beta_n(line, type);
                case THIRD:
                        return 0;
                default:
                        return core_l_get_beta_n(line, THIRD);
        }
}

typedef struct {
        double (*beta_n)(const apogee_rc_t *,
                         beta_ord_t);
        double (*alpha_n)(const apogee_rc_t *,
                          const unsigned int,
                          const double);
} coeff_t;

static coeff_t filler[TOTAL_QTY] = {
        { .beta_n = core_vr_get_beta_n,  .alpha_n = core_vr_get_alpha_n }, // VR_PART
        { .beta_n = united_l_get_beta_n, .alpha_n = core_l_get_alpha_n },  // L_PART
        { .beta_n = core_b_get_beta_n,   .alpha_n = core_b_get_alpha_n },  // B_PART
};

static void __fill_part_mnk_matrix(linear_equation_t *tmp,
                                   const apogee_rc_table_t *table,
                                   const unsigned int k)
{
        unsigned int j, i, n;
        unsigned int len = tmp->size;
        matrix_line_t *line = &g_matrix_line;

        memset(tmp->data, 0, sizeof(double) * len * len);
        memset(tmp->right, 0, sizeof(double) * len);

        for (j = 0; j < table->size; ++j) {
                for (i = 0; i < BETA_QTY + 1; ++i) {
                        line->_[i] = filler[k].beta_n(&table->data[j], i);
                }

                for (i = BETA_QTY + 1; i < len; ++i) {
                        line->_[i] = filler[k].alpha_n(&table->data[j],
                                                        i - BETA_QTY,
                                                        GET_TABLE_R0(table));
                }

                for (i = 0; i < len; ++i) {
                        double m = line->_[i];
                        for (n = 0; n < len; ++n) {
                                tmp->data[i * len + n] += line->_[n] * m / g_sq[k];
                        }

                        if (k == VR_PART) {
                                tmp->right[i] += table->data[j].v_helio * m / g_sq[k];
                                continue;
                        }

                        if (k == L_PART) {
                                tmp->right[i] += table->data[j].pm_l * m / g_sq[k];
                                continue;
                        }

                        if (k == B_PART) {
                                tmp->right[i] += table->data[j].pm_b * m / g_sq[k];
                                continue;
                        }
                }
        }
}

static void uni_fill_mnk_matrix(linear_equation_t *eq,
                                apogee_rc_table_t *table)
{
        unsigned int k;
        unsigned int len = eq->size; // A + thetas + omega_0 + v_sun + .. + w_sun

        memset(eq->data, 0, sizeof(double) * len * len);
        memset(eq->right, 0, sizeof(double) * len);

        linear_equation_t tmp;
        tmp.data = (double *)dv_alloc(sizeof(double) * len * len);
        tmp.right = (double *)dv_alloc(sizeof(double) * len);
        tmp.size = len;
        tmp.ord = eq->ord;

        for (k = 0; k < TOTAL_QTY; ++k) {
                __fill_part_mnk_matrix(&tmp, table, k);
                add_matrix_to_matrix(&tmp, eq);
        }
}

void precalc_errors_uni(apogee_rc_table_t *table,
                        const double limit)
{
        unsigned int i, j;
        double l_sq[TOTAL_QTY];

        for (j = 0; j < TOTAL_QTY; ++j) {
                l_sq[j] = sqrt(g_sq[j]);
        }

        for (i = 0; i < table->size; ++i) {
                for (j = 0; j < TOTAL_QTY; ++j) {
                        if (fabs(table->data[i].vsd[j]) / l_sq[j] > limit) {
                                table->data[i].pm_match = 0;
                        }
                }
        }
}


void precalc_errors_uni_sigma_0(apogee_rc_table_t *table,
				const double limit)
{
        unsigned int i, j;

        for (i = 0; i < table->size; ++i) {
                for (j = 0; j < TOTAL_QTY; ++j) {
			double s = sigma_for_k(j, &table->data[i], table->sigma_0);
                        if (fabs(table->data[i].vsd[j]) / s > limit) {
                                table->data[i].pm_match = 0;
                        }
                }
        }
}

double sigma_for_k(const unsigned int p,
                   const apogee_rc_t *line,
                   const double sigma_0)
{
        if (p == VR_PART) {
                return pow_double(sigma_0, 2) + pow_double(line->v_err, 2);
        }

        if (p == L_PART) {
                return pow_double(K_PM * sigma_0 / line->dist, 2) + pow_double(K_PM * line->pm_l_err, 2);
        }

        if (p == B_PART) {
                return pow_double(K_PM * sigma_0 / line->dist, 2) + pow_double(K_PM * line->pm_b_err, 2);
        }

        printf("%s: bad k %d\n", __func__, p);
        return 0;
}

static void uni_fill_mnk_matrix_nerr(linear_equation_t *eq,
                                     apogee_rc_table_t *table)
{
        unsigned int i, j, k, n;
        unsigned int len = eq->size; // A + thetas + omega_0 + v_sun + .. + w_sun
        matrix_line_t *line = &g_matrix_line;

        memset(eq->data, 0, sizeof(double) * len * len);
        memset(eq->right, 0, sizeof(double) * len);

        linear_equation_t tmp;
        tmp.data = (double *)dv_alloc(sizeof(double) * len * len);
        tmp.right = (double *)dv_alloc(sizeof(double) * len);
        tmp.size = len;
        tmp.ord = eq->ord;

        for (k = 0; k < TOTAL_QTY; ++k) {
                memset(tmp.data, 0, sizeof(double) * len * len);
                memset(tmp.right, 0, sizeof(double) * len);
                for (j = 0; j < table->size; ++j) {
                        double s = sigma_for_k(k, &table->data[j], table->sigma_0);
                        for (i = 0; i < BETA_QTY + 1; ++i) {
                                line->_[i] = filler[k].beta_n(&table->data[j], i);
                        }

                        for (i = BETA_QTY + 1; i < len; ++i) {
                                line->_[i] = filler[k].alpha_n(&table->data[j],
                                                                i - BETA_QTY,
                                                                GET_TABLE_R0(table));
                        }

                        for (i = 0; i < len; ++i) {
                                double m = line->_[i];
                                for (n = 0; n < len; ++n) {
                                        tmp.data[i * len + n] += line->_[n] * m / s;
                                }

                                if (k == VR_PART) {
                                        tmp.right[i] += table->data[j].v_helio * m / s;
                                }

                                if (k == L_PART) {
                                        tmp.right[i] += table->data[j].pm_l * m / s;
                                }

                                if (k == B_PART) {
                                        tmp.right[i] += table->data[j].pm_b * m / s;
                                }
                        }
                }

                add_matrix_to_matrix(&tmp, eq);
        }
}

double get_v_generic_from_uni(const linear_eq_solve_t *v,
                              const apogee_rc_t *line,
                              const double r_0,
                              const unsigned int type)
{
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY + 1; ++i) {
                mod_v += v->data[i] * filler[type].beta_n(line, i);
        }

        for (i = BETA_QTY + 1; i < v->size; ++i) {
                mod_v += v->data[i] * filler[type].alpha_n(line, i - BETA_QTY, r_0);
        }

        return mod_v;
}


static void _residuals_line(const linear_eq_solve_t *v,
                            apogee_rc_t *line,
                            const double r_0)
{
        unsigned int k;
        for (k = 0; k < TOTAL_QTY; ++k) {
                double mod_v = get_v_generic_from_uni(v, line, r_0, k);

                if (k == VR_PART) {
                        line->vsd[k] = line->v_helio - mod_v;
                        continue;
                }

                if (k == L_PART) {
                        line->vsd[k] = line->pm_l - mod_v;
                        continue;
                }

                if (k == B_PART) {
                        line->vsd[k] = line->pm_b - mod_v;
                        continue;
                }
        }
}

static double _chi_line(const linear_eq_solve_t *v,
                        apogee_rc_t *line,
                        const double r_0,
                        const double *sigma)
{
        unsigned int k;
        double res = 0;

        for (k = 0; k < TOTAL_QTY; ++k) {
                double mod_v = get_v_generic_from_uni(v, line, r_0, k);

                if (k == VR_PART) {
                        res += pow_double(line->v_helio - mod_v, 2) / g_sq[k];
                        continue;
                }

                if (k == L_PART) {
                        res += pow_double(line->pm_l - mod_v, 2) / g_sq[k];
                        continue;
                }

                if (k == B_PART) {
                        res += pow_double(line->pm_b - mod_v, 2) / g_sq[k];
                        continue;
                }
        }

        return res;
}



static double _residuals_line_nerr(const linear_eq_solve_t *v,
                                   apogee_rc_t *line,
                                   const double r_0,
                                   const double sigma_0)
{
        unsigned int i, k;
        double res = 0;
        line->eps = 0;

        for (k = 0; k < TOTAL_QTY; ++k) {
                double mod_v = 0;
                double s = sigma_for_k(k, line, sigma_0);
                for (i = 0; i < BETA_QTY + 1; ++i) {
                        mod_v += v->data[i] * filler[k].beta_n(line, i);
                }

                for (i = BETA_QTY + 1; i < v->size; ++i) {
                        mod_v += v->data[i] * filler[k].alpha_n(line, i - BETA_QTY, r_0);
                }

                if (k == VR_PART) {
                        res += pow_double(line->v_helio - mod_v, 2) / s;
                        line->vsd[k] = line->v_helio - mod_v;
                        continue;
                }

                if (k == L_PART) {
                        res += pow_double(line->pm_l - mod_v, 2) / s;
                        line->vsd[k] = line->pm_l - mod_v;
                        continue;
                }

                if (k == B_PART) {
                        res += pow_double(line->pm_b - mod_v, 2) / s;
                        line->vsd[k] = line->pm_b - mod_v;
                        continue;
                }
        }

        return res;
}

static double residuals_summary(const linear_eq_solve_t *solution,
                                apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i, j;
        const unsigned int n_free = table->size - solution->size - 1;
        for (i = 0; i < table->size; ++i) {
                _residuals_line(solution,
                                &table->data[i],
                                GET_TABLE_R0(table));
                for (j = 0; j < TOTAL_QTY; ++j)
                        table->sigma[j] += table->data[i].vsd[j] * table->data[i].vsd[j];
        }

        for (j = 0; j < TOTAL_QTY; ++j) {
                table->sigma[j] = table->sigma[j] / n_free;
        }

        for (i = 0; i < table->size; ++i) {
                sum += _chi_line(solution,
                                 &table->data[i],
                                 GET_TABLE_R0(table),
                                 table->sigma);
        }
        assert(sum > 0);
        return sum;
}

void precalc_vsd_to_dump(apogee_rc_table_t *table)
{
        unsigned int i, j;
	parser_t *cfg = get_parser();
        solution_mode_t mode = GET_SOLUTION_MODE(cfg);

        for (i = 0; i < table->size; ++i) {
		for (j = 0; j < TOTAL_QTY; ++j) {
			switch(mode) {
			case UNINAT_MODE:
			case FIND_SIGMA0_MODE:
				table->data[i].vsd[j] /= sigma_for_k(j,
								     &table->data[i],
								     cfg->sigma_0);
				break;
			case UNI_MODE:
				table->data[i].vsd[j] /= sqrt(g_sq[j]);
				break;
			default:
				// TODO: warnings
				return;
			}
		}
	}
}

static double residuals_summary_nerr(const linear_eq_solve_t *solution,
                                     apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                sum += _residuals_line_nerr(solution,
                                            &table->data[i],
                                            GET_TABLE_R0(table),
                                            table->sigma_0);
        }
        assert(sum > 0);
        return sum;
}



#if 0 // Here errors from LSE matrixes. Now it all from M.-K.
      // TODO: replace to deadcode
void uni_get_errors(opt_t *solution, apogee_rc_table_t *table)
{
        linear_equation_t m, invm;
        m.data = dv_alloc(sizeof(double) * solution->s.size * solution->s.size);
        m.right = dv_alloc(sizeof(double) * solution->s.size);
        m.size = solution->s.size;
        invm.data = dv_alloc(sizeof(double) * solution->s.size * solution->s.size);
        invm.size = solution->s.size;

        uni_fill_mnk_matrix(&m, table);
        inverse_and_diag(&m, &invm);

        solution->bounds = dv_alloc(sizeof(prec_t) * solution->s.size);
        unsigned int i;
        for (i = 0; i < solution->s.size; ++i) {
                solution->bounds[i].l =
                        get_error_mnk_estimated(invm.data[i * invm.size + i], invm.size + table->size + 1, solution->sq / (table->size + 1));
        }
}
#endif


static opt_t *united_solution(linear_equation_t *eq,
                              apogee_rc_table_t *table)
{
        opt_params_t params = {
                .residuals_summary = residuals_summary,
                .fill_mnk_matrix = uni_fill_mnk_matrix,
        };
        return opt_linear(eq, table, &params);
}

static opt_t *united_with_nerr_solution(linear_equation_t *eq,
                                        apogee_rc_table_t *table)
{
        opt_params_t params = {
                .residuals_summary = residuals_summary_nerr,
                .fill_mnk_matrix = uni_fill_mnk_matrix_nerr,
        };
        return opt_linear(eq, table, &params);
}



void uni_g_sd_init(const double *values)
{
        unsigned int i;
        for (i = 0; i < TOTAL_QTY; ++i) {
                g_sq[i] = values[i];
        }
}

opt_t *united_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        unsigned int dim = size + TOTAL_QTY + 1;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_t *opt = united_solution(&eq, table);

        update_table_R0(table, GET_SOLUTION_R0(opt));

        return opt;
}


void dump_united_sigma0_solution_profile(apogee_rc_table_t *table,
                                         unsigned int ord)
{
        unsigned int size = ord;
        unsigned int dim = size + TOTAL_QTY + 1;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_params_t params = {
                .residuals_summary = residuals_summary_nerr,
                .fill_mnk_matrix = uni_fill_mnk_matrix_nerr,
        };

        dump_profile(&eq, table, &params, "uni_profile.txt");
}

void dump_united_solution_profile(apogee_rc_table_t *table,
                                  unsigned int ord)
{
        unsigned int size = ord;
        unsigned int dim = size + TOTAL_QTY + 1;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_params_t params = {
                .residuals_summary = residuals_summary,
                .fill_mnk_matrix = uni_fill_mnk_matrix,
        };

        dump_profile(&eq, table, &params, "uni_profile.txt");
}

static void __find_r_0_bounds_f(apogee_rc_table_t *table,
			        opt_t *solution)
{
	parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        unsigned int dim = size + TOTAL_QTY + 1;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_params_t params = {
                .residuals_summary = residuals_summary_nerr,
                .fill_mnk_matrix = uni_fill_mnk_matrix_nerr,
        };

	find_r_0_bounds(table, solution, &params, &eq);
	dump_r0_bounds(solution);
}

static void __find_r_0_bounds_u(apogee_rc_table_t *table,
			        opt_t *solution)
{
	parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        unsigned int dim = size + TOTAL_QTY + 1;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_params_t params = {
                .residuals_summary = residuals_summary,
                .fill_mnk_matrix = uni_fill_mnk_matrix,
        };


	find_r_0_bounds(table, solution, &params, &eq);
	dump_r0_bounds(solution);
}

void dump_united_solution_r0_bounds(apogee_rc_table_t *table,
				    opt_t *solution)
{
	parser_t *cfg = get_parser();
        solution_mode_t mode = GET_SOLUTION_MODE(cfg);
	switch(mode) {
        case UNI_MODE:
		__find_r_0_bounds_u(table, solution);
		return;
        case FIND_SIGMA0_MODE:
		__find_r_0_bounds_f(table, solution);
		return;
        default:
		PR_WARN("Unknown united solution mode");
		return;
	}
}

opt_t *united_with_nature_errs_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        unsigned int size = cfg->ord;
        unsigned int dim = size + TOTAL_QTY + 1;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);
        table->sigma_0 = cfg->sigma_0;

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_t *opt = united_with_nerr_solution(&eq, table);

        update_table_R0(table, GET_SOLUTION_R0(opt));

        return opt;
}
