#include <math.h>
#include <stdio.h>

#include "math.h"
#include "asserts.h"
#include "mem.h"
#include "types.h"
#include "io.h"

#define HRD_JTRIG_CONST   9.606614

int __jtrigonometry(apogee_rc_table_t *a)
{
#ifdef DEBUG
        printf("%s: test distance = %lf\n",
                        __func__, get_R_distance(a->data, 8));
#endif
        return !(fabs(get_R_distance(a->data, 8) - HRD_JTRIG_CONST) < 1e-6);
}


int jtrigonometry()
{
        apogee_rc_table_t *table = read_table(INPUT_TABLE_FILE_NAME);
        return __jtrigonometry(table);
}

/**
 * From https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html
 */
int jmath()
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

        static double data_b[] = {
                1.0, 2.0, 3.0, 4.0
        };
        eq->right = data_b;

        linear_eq_solve_t *res = dv_alloc(sizeof(linear_eq_solve_t));
        res->size = 4;
        res->data = (double *)dv_alloc(sizeof(double) * res->size);

        solve(eq, res);

        unsigned int i;
        for (i = 0; i < size; ++i) {
                printf("%s: %f\n",
                                __func__, *(res->data + i));
        }

        printf("%s: kappa(5000) = %0.7lf\n",
                        __func__,
                        get_limit_by_eps(5000));
        assert(dv_factorial(0) == 1);
        assert(dv_factorial(1) == 1);
        assert(dv_factorial(2) == 2);
        assert(dv_factorial(3) == 6);
        assert(dv_factorial(4) == 24);
        return 0;
}
