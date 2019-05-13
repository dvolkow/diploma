#ifndef CORE_B_H
#define CORE_B_H  1
#include "types.h"
#include "opt.h"

double core_b_get_beta_n(const apogee_rc_t *, beta_ord_t);
double core_b_get_alpha_n(const apogee_rc_t *,
                          const unsigned int,
                          const double);
opt_t *core_b_entry(apogee_rc_table_t *);
void precalc_errors_mu_b(apogee_rc_table_t *,
                         const double);
void get_partial_b_solution(apogee_rc_table_t *table);
#endif // CORE_B_H

