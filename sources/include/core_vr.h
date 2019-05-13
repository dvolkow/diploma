#ifndef CORE_VR_H
#define CORE_VR_H  1
#include "types.h"
#include "opt.h"

#define BETA_QTY_FIX    2


opt_t *core_vr_entry(apogee_rc_table_t *);
double __core_vr_get_beta_n(const apogee_rc_t *line,
                            beta_ord_t type);
void precalc_errors_vr(apogee_rc_table_t *,
                       const double);
#endif // CORE_VR_H
