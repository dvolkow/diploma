#ifndef TRIGONOMETRY_H
#define TRIGONOMETRY_H

#include "types.h"

/*
 * By XEphem and Elwood Downey (J2000.0):
 */
#define L_GAL                   33
#define ALPHA_N                 193
#define DELTA_N                 27.17

#define K_PM                    4.7406


double mu_b_from_pa_dec_pm(const apogee_rc_t *);
double mu_l_from_pa_dec_pm(const apogee_rc_t *);
double errors_ecliptic_to_gal(const apogee_rc_t *);

#endif  // TRIGONOMETRY_H
