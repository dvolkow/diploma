#include "asserts.h"
#include "math.h"
#include "types.h"
#include "mem.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#ifdef DEBUG
#include "debug.h"
#endif


#define to_gsl_matrix(a)        \
        gsl_matrix_view_array((a)->data, (a)->size, (a)->size)

#define to_gsl_vector(a)        \
        gsl_vector_view_array((a)->right, (a)->size)

static double __factorial_storage[PRECACHED_FACTORIAL_LEN] = {
        1, 1, 2, 6, 24, 
        120, 720, 5040,
        40320, 362880
};



void *make_linear_struct(double *data, int size, 
                                linear_type_t type)
{
        void *block;

        switch(type) {
        case EQUATION:
                block = dv_alloc(sizeof(linear_equation_t)); 
                break;
        case SOLUTION:
                block = dv_alloc(sizeof(linear_eq_solve_t)); 
                break;
        default:
                printf("%s: warning: unknown type\n",
                                __FUNCTION__);
                block = NULL;
                break;
        }
        return block;
}

void gsl_vector_copy_to(linear_eq_solve_t *x, const gsl_vector *src)
{
        unsigned int i;
        for (i = 0; i < x->size; ++i) {
                x->data[i] = gsl_vector_get(src, i);
        }
}

void solve(linear_equation_t *eq, linear_eq_solve_t *x)
{
        gsl_vector *gx = gsl_vector_alloc(eq->size);
        int s;

        gsl_matrix_view m = to_gsl_matrix(eq);
        gsl_vector_view gb = to_gsl_vector(eq);
        gsl_permutation *p = gsl_permutation_alloc(eq->size);

        gsl_linalg_LU_decomp(&m.matrix, p, &s);
        gsl_linalg_LU_solve(&m.matrix, p, &gb.vector, gx);

        gsl_vector_copy_to(x, gx);
        gsl_permutation_free(p);
        gsl_vector_free(gx);
}

int dv_factorial(const int n)
{
        return __factorial_storage[n];
}

double dot_prod(double *a, double *b, int size) 
{
        int i = 0;
        double res = 0;
        for (i = 0; i < size; ++i) {
                res += a[i] * b[i];
        }
        return res; 
}
