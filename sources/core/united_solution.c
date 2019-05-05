#include <math.h>
#include <assert.h>

#include "types.h"
#include "core.h"
#include "core_l.h"
#include "core_b.h"
#include "unicore.h"

static matrix_line_t g_matrix_line;
static double g_sq[TOTAL_QTY];


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

static void uni_fill_mnk_matrix(linear_equation_t *eq,
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
                        for (i = 0; i < BETA_QTY + 1; ++i) {
                                line->_[i] = filler[k].beta_n(&table->data[j], i);
                        }

                        for (i = BETA_QTY + 1; i < len; ++i) {
                                line->_[i] = filler[k].alpha_n(&table->data[j], i - BETA_QTY, table->r_0);
                        }

                        for (i = 0; i < len; ++i) {
                                double m = line->_[i];
                                for (n = 0; n < len; ++n) {
                                        tmp.data[i * len + n] += line->_[n] * m / g_sq[k];
                                }

                                if (k == VR_PART) {
                                        tmp.right[i] += table->data[j].v_helio * m / g_sq[k];
                                }

                                if (k == L_PART) {
                                        tmp.right[i] += table->data[j].pm_l * m / g_sq[k];
                                }

                                if (k == B_PART) {
                                        tmp.right[i] += table->data[j].pm_b * m / g_sq[k];
                                }
                        }
                }

                add_matrix_to_matrix(&tmp, eq);
        }
}

static double _residuals_line(const linear_eq_solve_t *v,
                             apogee_rc_t *line,
                             const double r_0)
{
        unsigned int i, k;
        double res = 0;
        for (k = 0; k < TOTAL_QTY; ++k) {
                double mod_v = 0;
                for (i = 0; i < BETA_QTY + 1; ++i) {
                        mod_v += v->data[i] * filler[k].beta_n(line, i);
                }

                for (i = BETA_QTY + 1; i < v->size; ++i) {
                        mod_v += v->data[i] * filler[k].alpha_n(line, i - BETA_QTY, r_0);
                }

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

        //line->eps = fabs(line->v_helio - mod_v);
        return res;
}


static double residuals_summary(const linear_eq_solve_t *solution, 
                                const apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                sum += _residuals_line(solution, &table->data[i],
                                        table->r_0);
        }
        assert(sum > 0);
        return sum;
}


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



static opt_t *united_solution(linear_equation_t *eq,
                              apogee_rc_table_t *table)
{
        opt_params_t params = {
                .residuals_summary = residuals_summary,
                .fill_mnk_matrix = uni_fill_mnk_matrix,
        };
        return opt_linear(eq, table, &params);
}


void uni_g_sd_init(const double *values)
{
        unsigned int i;
        for(i = 0; i < TOTAL_QTY; ++i) {
                g_sq[i] = values[i];
        }
}

opt_t *united_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        int size = cfg->ord;
        int dim = size + TOTAL_QTY + 1;
        double *matrix = (double *)dv_alloc(sizeof(double) * dim * dim);

#ifdef DEBUG
        printf("SD^2:\n");
        printf("r_0: %lf\n", g_sq[VR_PART]);
        printf("l: %lf\n", g_sq[L_PART]);
        printf("b: %lf\n", g_sq[B_PART]);
#endif

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * dim),
                .size = dim,
                .ord = size
        };

        opt_t *opt = united_solution(&eq, table);
        uni_get_errors(opt, table);
        table->r_0 = opt->r_0;

        return opt;
}
