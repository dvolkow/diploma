#ifndef GENERATORS_H
#define GENERATORS_H    1

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>

#include "opt.h"
#include "types.h"


#define PREALLOC_TABLE_SIZE     32768

typedef enum {
        MT
      , ANOTHER
} gen_mode_t;

typedef struct {
        opt_t *(*f_entry)(apogee_rc_table_t *);
        double (*f_point_by_solution)(const opt_t *,
                                      const double);
        void (*f_table_by_solution)(const opt_t *,
                                    const apogee_rc_table_t *,
                                    apogee_rc_table_t *);
        unsigned int count;
	char *mul_unfres_name;
} mk_params_t;

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
        return (unsigned long)(tv.tv_sec) + (unsigned long)(tv.tv_usec);
}

opt_t *monte_carlo_entry(const opt_t *,
                         const apogee_rc_table_t *,
                         const mk_params_t *);

void fill_table_by_uni_solution(const opt_t *,
                                const apogee_rc_table_t *,
                                apogee_rc_table_t *);
void fill_table_by_uni_solution_sigma_0(const opt_t *,
					const apogee_rc_table_t *,
					apogee_rc_table_t *);

void fill_table_by_vr_solution(const opt_t *,
                               const apogee_rc_table_t *,
                               apogee_rc_table_t *);
void fill_table_by_b_solution(const opt_t *,
                               const apogee_rc_table_t *,
                               apogee_rc_table_t *);
void fill_table_by_l_solution(const opt_t *,
                               const apogee_rc_table_t *,
                               apogee_rc_table_t *);
double *gen_vector_by_bounds_uni(const gsl_rng *r,
                                 const double l,
                                 const double h,
                                 const unsigned int size);

void generate(void);
apogee_rc_table_t *create_apogee_rc_table_by_size(unsigned int);
apogee_rc_table_t *generic_table(void);

int random_seed_init(void);
void random_seed_exit(void);
#endif // GENERATORS_H
