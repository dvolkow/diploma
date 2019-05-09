#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "math.h"
#include "mem.h"
#include "types.h"
#include "core.h"
#include "opt.h"
#include "unicore.h"

#define BETA_QTY_FIX    2

static double w_0 = 0;
static matrix_line_t g_matrix_line;

static double __core_vr_get_beta_n(const apogee_rc_t *line,
                                 beta_ord_t type)
{
        switch (type) {
                case FIRST:
                        return -cos(line->l) * cos(line->b);
                case SECOND:
                        return -sin(line->l) * cos(line->b);
                case THIRD:
                        // fixed parameter:
                        return w_0 * (-sin(line->b));
                default:
#ifdef DEBUG
                        printf("%s: type error!\n", __func__);
#endif
                        return 0;
        }
}

static void core_vr_fill_mnk_matrix(linear_equation_t *eq,
                                    apogee_rc_table_t *table)
{
        unsigned int i, j, k;
        unsigned int len = eq->size;
        matrix_line_t *line = &g_matrix_line;
#ifdef DEBUG
        assert(len <= BETA_QTY_FIX + MAX_ORDER_SOLUTION);
#endif

        memset(eq->data, 0, sizeof(double) * eq->size * eq->size);
        memset(eq->right, 0, sizeof(double) * eq->size);

        for (j = 0; j < table->size; ++j) {
                for (i = 0; i < BETA_QTY_FIX; ++i) {
                        line->_[i] = __core_vr_get_beta_n(&table->data[j], i);
                }

                for (i = BETA_QTY_FIX; i < eq->size; ++i) {
                        line->_[i] = core_vr_get_alpha_n(&table->data[j], i - BETA_QTY_FIX + 1, table->r_0);
                }

                for (i = 0; i < len; ++i) {
                        double m = line->_[i];
                        for (k = 0; k < len; ++k) {
                                eq->data[i * len + k] += line->_[k] * m;
                        }
                        eq->right[i] += (table->data[j].v_helio - __core_vr_get_beta_n(&table->data[j], THIRD)) * m;
                }
        }
}

static double core_vr_get_mod_v(const opt_t *solution,
                                const apogee_rc_t *line)
{
        const double r_0 = solution->r_0;
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY_FIX; ++i) {
                mod_v += solution->s.data[i] * __core_vr_get_beta_n(line, i);
        }

        for (i = BETA_QTY_FIX; i < solution->s.size; ++i) {
                mod_v += core_vr_get_alpha_n(line, i - BETA_QTY_FIX + 1, r_0) *
                                solution->s.data[i];
        }

        // add fixed parameter:
        return mod_v + __core_vr_get_beta_n(line, THIRD);
}

static double _residuals_line(const linear_eq_solve_t *v,
                             apogee_rc_t *line,
                             const double r_0)
{
        double mod_v = 0;
        unsigned int i;
        for (i = 0; i < BETA_QTY_FIX; ++i) {
                mod_v += v->data[i] * __core_vr_get_beta_n(line, i);
        }
        for (i = BETA_QTY_FIX; i < v->size; ++i) {
                mod_v += core_vr_get_alpha_n(line, i - BETA_QTY_FIX + 1, r_0) * v->data[i];
        }

        line->eps = fabs(line->v_helio - mod_v - __core_vr_get_beta_n(line, THIRD));
        return pow_double(line->v_helio - mod_v - __core_vr_get_beta_n(line, THIRD), 2);
}


static double residuals_summary(const linear_eq_solve_t *solution, 
                                apogee_rc_table_t *table)
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

void precalc_errors_vr(apogee_rc_table_t *table,
                       const double limit)
{
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                if (table->data[i].eps / table->sigma[VR_PART] > limit) {
                        table->data[i].pm_match = 0;
                }
        }
}

static opt_t *core_vr_get_solution(linear_equation_t *eq,
                                   apogee_rc_table_t *table)
{
        opt_params_t params = {
                .residuals_summary = residuals_summary,
                .fill_mnk_matrix = core_vr_fill_mnk_matrix
        };
        return opt_linear(eq, table, &params);
}


opt_t *core_vr_entry(apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        int size = cfg->ord;
        double *matrix = (double *)dv_alloc(sizeof(double) * (size + BETA_QTY_FIX) *
                                                   (size + BETA_QTY_FIX));

        // fixed parameter:
        w_0 = table->w_sun;

        linear_equation_t eq = {
                .data = matrix,
                .right = (double *)dv_alloc(sizeof(double) * (size + BETA_QTY_FIX)),
                .size = size + BETA_QTY_FIX,
                .ord = size
        };

        opt_t *opt = core_vr_get_solution(&eq, table);
        opt->sq = sqrt(opt->sq / (table->size - eq.size - 1));
        table->r_0 = opt->r_0;
        table->sigma[VR_PART] = opt->sq;

        dump_core_vr_solution(opt);
        return opt;
}
