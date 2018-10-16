#ifndef GENERATORS_H
#define GENERATORS_H    1

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>

typedef enum {
        MT
      , ANOTHER
} gen_mode_t;

double *gen_vector_by_mean_and_sd(const gsl_rng *r, 
                                const double mean, 
                                 const double sigma,
                                 const unsigned int size);

gsl_rng *dv_rand_acquire(gen_mode_t mode);
void dv_rand_release(gsl_rng *r);

static inline unsigned long int dv_random_seed(void)
{
        struct timeval tv;
        gettimeofday(&tv, 0);
        return (tv.tv_sec + tv.tv_usec);
}

double *gen_vector_by_bounds_uni(const gsl_rng *r,
                                 const double l,
                                 const double h,
                                 const unsigned int size);

void generate(void);

int random_seed_init(void);
void random_seed_exit(void);
#endif // GENERATORS_H
