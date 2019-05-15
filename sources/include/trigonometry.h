#ifndef TRIGONOMETRY_H
#define TRIGONOMETRY_H

#include "types.h"

/*
 * By XEphem and Elwood Downey (J2000.0):
 */
#define L_GAL                   32.93192
#define ALPHA_N                 192.85948
#define DELTA_N                 27.12825

#define K_PM                    4.7406


double errors_ecliptic_to_gal_mu_l(const apogee_rc_t *);
double errors_ecliptic_to_gal_mu_b(const apogee_rc_t *);

double mu_b_from_pa_dec_pm_II(const apogee_rc_t *);
double mu_l_from_pa_dec_pm_II(const apogee_rc_t *);
#endif  // TRIGONOMETRY_H
