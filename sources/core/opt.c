#include <math.h>
#include "asserts.h"

#include "types.h"
#include "math.h"



static double residuals_line(const linear_equation_t *eq,
                             const linear_eq_solve_t *v,
                             const apogee_rc_t *line)
{
        double mod_v = dot_prod(eq->data, v->data, v->size);
        return (line->v_helio - mod_v) * (line->v_helio - mod_v);
}

static double residuals_summary(const linear_equation_t *eq, 
                                const linear_eq_solve_t *v, 
                                const apogee_rc_table_t *table)
{
        double sum = 0;
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                sum += residuals_line(eq, v, &table->data[i]);
        }
        return sum;
}

static void fill_right_vector(const apogee_rc_table_t *table,
                                linear_eq_solve_t *r)
{
        unsigned int i;
        assert(table->size == r->size);
        for (i = 0; i < table->size; ++i) {
                r->data[i] = table->data[i].v_helio;
        }
}

typedef struct {
        linear_eq_solve_t s;
        double r_0;
        double sq;
} opt_t;

linear_eq_solve_t *opt_linear(linear_equation_t *eq, 
                                apogee_rc_table_t *table) 
{
        unsigned int i;
        linear_eq_solve_t *s = dv_alloc(sizeof(linear_eq_solve_t));
        s->data = dv_alloc(sizeof(double) * eq->size);

        linear_eq_solve_t *r = dv_alloc(sizeof(linear_eq_solve_t));
        r->data = dv_alloc(sizeof(double) * table->size);
        r->size = table->size;
        fill_right_vector(table, r);

        double low_r = LOWER_BOUND_R0;
        double high_r = UPPER_BOUND_R0;
        double step = (high_r - low_r) / 32;

        table->r_0 = low_r;

        fill_mnk_matrix_vr(eq->data, eq->ord, table);
        solve(eq, r, s); 
        double sq = residuals_summary(eq, s, table);

        opt_t opt_params = {
                .s = *s,
                .r_0 = table->r_0,
                .sq = sq
        };

        while (step > SEARCH_PRECISION) {
                while (table->r_0 < high_r) {
                        table->r_0 += step;
                        fill_mnk_matrix_vr(eq->data, eq->ord, table);
                        solve(eq, r, s);
                        double sq_tmp = residuals_summary(eq, s, table);
                        if (sq_tmp < sq) {
                                opt_params.s = *s;
                                opt_params.r_0 = table->r_0;
                                opt_params.sq = sq_tmp;
                        }
                }

                low_r = opt_params.r_0 - step;
                high_r = opt_params.r_0 + step;
                step /= STEP_DIVISOR;
        }

        return s;
}

