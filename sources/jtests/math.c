#include "math.h"
#include "asserts.h"
#include "mem.h"
#include "types.h"

void jtrigonometry(apogee_rc_t *a)
{
        printf("%s: test distance = %lf\n",
                        __FUNCTION__, get_R_distance(a, 8));
}

/**
 * From https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html
 */
void jmath()
{
        int size = 4;

        linear_equation_t *eq = dv_alloc(sizeof(linear_equation_t));
        eq->size = 4;
        static double data_eq[] = {
                0.18, 0.60, 0.57, 0.96,
                0.41, 0.24, 0.99, 0.58,
                0.14, 0.30, 0.97, 0.66,
                0.51, 0.13, 0.19, 0.85
        };
        eq->data = data_eq;

        linear_eq_solve_t *b = dv_alloc(sizeof(linear_eq_solve_t));
        b->size = 4;
        static double data_b[] = {
                1.0, 2.0, 3.0, 4.0
        };
        b->data = data_b;

        linear_eq_solve_t *res = dv_alloc(sizeof(linear_eq_solve_t));
        res->size = 4;
        res->data = (double *)dv_alloc(sizeof(double) * res->size);

        solve(eq, b, res);

        unsigned int i;
        for (i = 0; i < size; ++i) {
                printf("%s: %f\n",
                                __FUNCTION__, *(res->data + i));
        }

}
