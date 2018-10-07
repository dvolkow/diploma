#include <math.h>

#include "math.h"
#include "asserts.h"
#include "mem.h"
#include "types.h"
#include "io.h"

#define HRD_JTRIG_CONST  7.821446

int __jtrigonometry(apogee_rc_table_t *a)
{
        /*
        printf("%s: test distance = %lf\n",
                        __FUNCTION__, get_R_distance(a->data, 8));
        */
        return !(fabs(get_R_distance(a->data, 8) - HRD_JTRIG_CONST) < 1e-7);
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

        return 0;
}
