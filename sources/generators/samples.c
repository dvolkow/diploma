#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "math.h"
#include "mem.h"
#include "generators.h"

static gsl_rng *g_gauss_rg_p;


double *gen_vector_by_mean_and_sd(const gsl_rng *r, 
                                  const double mean, 
                                  const double sigma,
                                  const unsigned int size)
{
        unsigned int i;
        double *res = dv_alloc(size);
        for (i = 0; i < size; ++i) {
                res[i] = mean + gsl_ran_gaussian(r, sigma);
        }

        return res;
}

static gsl_rng *__dv_rand_acquire(unsigned long int seed)
{
        gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(r, seed);

        return r;
}

gsl_rng *dv_rand_acquire(void)
{
#define get_random_double_by_memory()  \
        (*(double *)dv_mm_get_current_top())
        printf("%s: get_random_double_by_memory = %lf\n",
                        __func__, get_random_double_by_memory());
        return __dv_rand_acquire((unsigned long int)(gsl_rng_uniform(g_gauss_rg_p) 
                                        * 1000 + get_random_double_by_memory()));
}

void dv_rand_release(gsl_rng *r)
{
        gsl_rng_free(r);
}


int random_seed_init()
{
        g_gauss_rg_p = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(g_gauss_rg_p, dv_random_seed());
        return 0;
}

void random_seed_exit()
{
        dv_rand_release(g_gauss_rg_p);
}
