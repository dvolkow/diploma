#ifndef OPT_H
#define OPT_H   1

#include "math.h"
#include "types.h"

#define BOUNDS_R0_PREC	(1e-3)

typedef enum {
        LOWER,
        UPPER
} bound_t;


typedef struct {
        double (*residuals_summary)(const linear_eq_solve_t *,
                                    apogee_rc_table_t *);
        void (*fill_mnk_matrix)(linear_equation_t *,
                                apogee_rc_table_t *);
} opt_params_t;

opt_t *opt_linear(linear_equation_t *,
                  apogee_rc_table_t *,
                  opt_params_t *);

double lower_bound_search(linear_equation_t *eq,
                                apogee_rc_table_t *table,
                                double r_0);
double upper_bound_search(linear_equation_t *eq,
                                apogee_rc_table_t *table,
                                double r_0);
void get_errors(opt_t *, apogee_rc_table_t *table);
double opt_residuals_summary(const linear_eq_solve_t *,
                             apogee_rc_table_t *);

double get_beta_n(const apogee_rc_t *line, beta_ord_t type);
double get_mod_vr(const opt_t *solution,
                  const apogee_rc_t *line);
double get_mod_b(const opt_t *solution,
                 const apogee_rc_t *line);
double get_mod_l(const opt_t *solution,
                 const apogee_rc_t *line);


opt_t *exception_algorithm(apogee_rc_table_t *,
                           opt_t *(*f)(apogee_rc_table_t *),
                           void (*precalc_errors)(apogee_rc_table_t *,
                                                  const double));


void find_r_0_bounds(apogee_rc_table_t *,
		     opt_t *,
		     opt_params_t *,
		     linear_equation_t *);
#endif // OPT_H
