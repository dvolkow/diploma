#include <math.h>
#include "asserts.h"

#include "types.h"
#include "math.h"
#ifdef DEBUG
#include "debug.h"
#endif // DEBUG


double get_beta_n(const apogee_rc_t *line, beta_ord_t type);
double get_alpha_n(const apogee_rc_t *line, 
                   const double r_0,
                   const int n);

void fill_mnk_matrix_vr(linear_equation_t *eq,
                         apogee_rc_table_t *table);
/**
 * TODO: generic
 */
static double residuals_line(const linear_equation_t *eq,
                             const linear_eq_solve_t *v,
                             const apogee_rc_t *line,
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

opt_t *opt_linear(linear_equation_t *eq, 
                                apogee_rc_table_t *table) 
{
        unsigned int i;
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
                .s = s,
                .r_0 = table->r_0,
                .sq = sq
        };

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
                                        __FUNCTION__, sq_tmp, table->r_0);
        #endif
                        if (sq_tmp < opt_params.sq) {
                                opt_params.s = s;
                                opt_params.r_0 = table->r_0;
                                opt_params.sq = sq_tmp;
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

