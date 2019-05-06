#include "asserts.h"
#include "math.h"
#include "types.h"
#include "mem.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_erf.h>

#ifdef DEBUG
#include "debug.h"
#endif


#define to_gsl_matrix(a)        \
        gsl_matrix_view_array((a)->data, (a)->size, (a)->size)

#define to_gsl_vector(a)        \
        gsl_vector_view_array((a)->right, (a)->size)

static double __factorial_storage[PRECACHED_FACTORIAL_LEN];

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
                                __func__);
                block = NULL;
                break;
        }
        return block;
}

static void gsl_vector_copy_to(linear_eq_solve_t *x, const gsl_vector *src)
{
        unsigned int i;
        for (i = 0; i < x->size; ++i) {
                x->data[i] = gsl_vector_get(src, i);
        }
}


/*
 * Ax = B linear systems solver
 */
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


static void copy_diagonal(linear_equation_t *dst, 
                                const gsl_matrix *src)
{
        unsigned int i;
        for (i = 0; i < dst->size; ++i) {
                dst->data[i * dst->size + i] =
                        gsl_matrix_get(src, i, i);
        }
}

void add_matrix_to_matrix(const linear_equation_t *src,
                          linear_equation_t *dst)
{
        unsigned int i, j;
        assert(src->size == dst->size);
        const unsigned int len = dst->size;

        for (i = 0; i < len; ++i) {
                for (j = 0; j < len; ++j) {
                        dst->data[i * len + j] +=
                                src->data[i * len + j];
                }

                dst->right[i] += src->right[i];
        }
}

/*
 * Get inverse matrix for @eq->data, and get only
 * diagonal elems. 
 */
void inverse_and_diag(linear_equation_t *eq, linear_equation_t *res)
{
        gsl_matrix_view m = to_gsl_matrix(eq);
        gsl_matrix_view invm = to_gsl_matrix(res);
        gsl_permutation *p = gsl_permutation_alloc(eq->size);

        int s;
        gsl_linalg_LU_decomp(&m.matrix, p, &s);
        gsl_linalg_LU_invert(&m.matrix, p, &invm.matrix);

        copy_diagonal(res, &invm.matrix);

        gsl_permutation_free(p);
}


/*
 * Use precalculated values
 */
double dv_factorial(const int n)
{
        assert(n < PRECACHED_FACTORIAL_LEN);
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

double get_error_mnk_estimated(const double p, __attribute__((__unused__)) const int nfree,
                                const double sd)
{
        return sqrt(p * sd);
}


/* Statistics: */
double get_median(const double *data, const size_t size)
{
        return gsl_stats_median_from_sorted_data(data, 1, size);
}

double get_mean(const double *data, const size_t size)
{
        return gsl_stats_mean(data, 1, size); 
}

double get_sd(const double *data, const size_t size)
{
        return gsl_stats_sd(data, 1, size);
}

/* Error functions: */
static double __psi(const double kappa)
{
        return 1 - 2 * gsl_sf_erf_Q(kappa);
}

double get_limit_by_eps(const unsigned int size)
{
        const double left = 1 -  1. / size; 
        const double step = 1e-4;
        double kappa = 2;
        while (__psi(kappa) < left) {
                kappa += step;
        }
        return kappa - step;
}

int math_init(void) 
{
        __factorial_storage[0] = 1;
        int i;
        for (i = 1; i < PRECACHED_FACTORIAL_LEN; ++i) {
                __factorial_storage[i] = __factorial_storage[i - 1] * i;
        }

        return 0;
}

void math_exit() {}
