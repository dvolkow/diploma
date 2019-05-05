#ifndef CORE_B_H  
#define CORE_B_H  1
#include "types.h"
#include "opt.h"

double core_b_get_beta_n(const apogee_rc_t *, beta_ord_t);
double core_b_get_alpha_n(const apogee_rc_t *,
                          const unsigned int,
                          const double);
opt_t *core_b_entry(apogee_rc_table_t *);
#endif // CORE_B_H

