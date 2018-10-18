#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <assert.h>

#include "math.h"
#include "mem.h"
#include "opt.h"
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

static inline double __gen_gauss_d(const gsl_rng *r, 
                                   const double mean, 
                                   const double sigma)
{
        return mean + gsl_ran_gaussian(r, sigma);
}

static inline double __gen_flat_d(const gsl_rng *r,
                                 const double l,
                                 const double h)
{
        return gsl_ran_flat(r, l, h);
}

static double __gen_dist_d(const gsl_rng *r,
                           const double mean, 
                           const double sd)
{
        double res;
        do {
                res = fabs(__gen_gauss_d(r, mean, sd));
        } 
        while (res > 12);
        return res;
}

static double __gen_b_d(const gsl_rng *r,
                           const double mean, 
                           const double sd)
{
#define CUT_BY_B        80
        double res; 
        do {
                res = __gen_gauss_d(r, mean, sd);
        } 
        while (fabs(res) >= CUT_BY_B);

        return res;
}

apogee_rc_table_t *gen_table_by_solution(const opt_t *solution)
{
        apogee_rc_table_t *table = dv_alloc(sizeof(apogee_rc_table_t));
        gsl_rng *rng = dv_rand_acquire(MT);
        assert(rng != NULL);


        const dsize_t ssize = solution->size;
        table->r_0 = solution->r_0;
        table->size = ssize;
        table->data = dv_alloc(sizeof(apogee_rc_t) * ssize);

        dsize_t i;
#define LLOW    0
#define LHIGH   360

#define BLOW    -90
#define BHIGH   90

#define RLOW    0.39
#define RHIGH   12.83    

        for (i = 0; i < ssize; ++i) {
                table->data[i].id = i;
                table->data[i].l = deg_to_rad(__gen_flat_d(rng, LLOW, LHIGH));
                table->data[i].b = deg_to_rad(__gen_b_d(rng, 0, 15)); //__gen_flat_d(rng, BLOW, BHIGH);
                table->data[i].dist = __gen_dist_d(rng, 3, 1.5); //__gen_flat_d(rng, RLOW, RHIGH);
                table->data[i].v_helio = __gen_gauss_d(rng, get_mod_vr(solution, 
                                                                &table->data[i]), solution->sq);
        }

        dv_rand_release(rng);
        return table;
}

double *gen_vector_by_bounds_uni(const gsl_rng *r,
                                 const double l,
                                 const double h,
                                 const unsigned int size)
{
        unsigned int i;
        double *res = dv_alloc(size);

        for (i = 0; i < size; ++i) {
                res[i] = __gen_flat_d(r, l, h);
        }

        return res;
}

static gsl_rng *__dv_rand_acquire(const unsigned long int seed,
                                  gen_mode_t mode)
{
        gsl_rng *r;
        switch (mode) {
                case MT:
                default:
                        r = gsl_rng_alloc(gsl_rng_mt19937);
                        break;
        }
        gsl_rng_set(r, seed);

        return r;
}

gsl_rng *dv_rand_acquire(gen_mode_t mode)
{
#define get_random_double_by_memory()  \
        (*(double *)dv_mm_get_current_top())
#ifdef DEBUG
        printf("%s: get_random_double_by_memory = %lf\n",
                        __func__, get_random_double_by_memory());
#endif
        return __dv_rand_acquire((unsigned long int)(gsl_rng_uniform(g_gauss_rg_p) 
                                        * 1000 + get_random_double_by_memory()), mode);
}

void dv_rand_release(gsl_rng *r)
{
        gsl_rng_free(r);
}


void generate(void)
{
        parser_t *cfg = get_parser();

        opt_t *solution = read_solution(cfg->input_file_name);
        if (solution == NULL) {
                printf("%s: failure when read solution\n", __func__);
                return;
        }

        apogee_rc_table_t *table = gen_table_by_solution(solution);
        dump_table(table);
        dump_objects_xyz(table, table->size);
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
