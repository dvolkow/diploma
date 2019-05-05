#ifndef CORE_L_H  
#define CORE_L_H  1
#include "types.h"

double core_l_get_beta_n(const apogee_rc_t *, beta_ord_t);
double core_l_get_alpha_n(const apogee_rc_t *,
                          const unsigned int,
                          const double);

opt_t *core_l_entry(apogee_rc_table_t *);
#endif // CORE_L_H
