#include <math.h>

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

linear_eq_solve_t *opt_linear(linear_equation_t *eq) 
{
        unsigned int i;
        linear_eq_solve_t *s = dv_alloc(sizeof(linear_eq_solve_t));
        return s;
}

