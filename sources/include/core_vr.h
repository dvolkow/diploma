#ifndef CORE_VR_H  
#define CORE_VR_H  1
#include "types.h"
#include "opt.h"


opt_t *core_vr_entry(apogee_rc_table_t *);
void precalc_errors_vr(apogee_rc_table_t *,
                       const double);
#endif // CORE_VR_H
