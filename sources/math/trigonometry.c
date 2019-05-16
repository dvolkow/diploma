#include <math.h>
#include <assert.h>

#include "types.h"
#include "math.h"
#include "mem.h"
#include "trigonometry.h"

#ifdef DEBUG
#include "debug.h"
#endif


double get_R_distance(const apogee_rc_t *line, const double r_0)
{
        double dist = line->dist;
#ifdef DEBUG
        assert(dist > 0);
        assert(r_0 > 0);
#endif
        double root = dist * dist * line->cos_b * line->cos_b +
                r_0 * r_0 - 2 * r_0 * dist * line->cos_l * line->cos_b;

#ifdef DEBUG
        assert(root >= 0);
#endif
        return sqrt(root);
}


static double get_x(const apogee_rc_t *line)
{
        return line->dist * line->cos_l * line->cos_b;
}

static double get_y(const apogee_rc_t *line)
{
        return line->dist * line->sin_l * line->cos_b;
}

static double get_z(const apogee_rc_t *line)
{
        return line->dist * line->sin_b;
}

point_t *get_point(const apogee_rc_t *line)
{
        point_t *p = dv_alloc(sizeof(point_t));
        p->x = get_x(line);
        p->y = get_y(line);
        p->z = get_z(line);

        return p;
}

static double get_sin_phi(const apogee_rc_t *line)
{
        return cos(deg_to_rad(DELTA_N)) * sin(line->ra - deg_to_rad(ALPHA_N)) / cos(line->b);
}

static double get_cos_phi(const apogee_rc_t *line)
{
        return (sin(deg_to_rad(DELTA_N)) - sin(line->b) * sin(line->dec)) / (cos(line->b) * cos(line->dec));
}

static double get_sin_b(const apogee_rc_t *line)
{
        return sin(line->dec) * cos(deg_to_rad(90 - DELTA_N)) -
                cos(line->dec) * sin(line->ra - deg_to_rad(ALPHA_N + 90)) * sin(deg_to_rad(90 - DELTA_N));
}

static double get_sin_phi_II(const apogee_rc_t *line)
{
        return (cos(line->dec) * sin(line->ra - deg_to_rad(ALPHA_N + 90)) * cos(deg_to_rad(90 - DELTA_N)) +
                                sin(line->dec) * sin(deg_to_rad(90 - DELTA_N))) / cos(asin(get_sin_b(line)));
}

static double get_cos_phi_II(const apogee_rc_t *line)
{
        return cos(line->dec) * cos(line->ra - deg_to_rad(ALPHA_N + 90)) / cos(asin(get_sin_b(line)));
}


#if 0 // Formulae need correction. Use numerical translation.
      // aka mu_l_from_pa_dec_pm_II etc.
/*
 * Conversion between galactical and ecliptical velocities:
 *
 * 1. Dotting for equations for angles
 * 2. Substitution and conversion coordinates
 */
double mu_l_from_pa_dec_pm(const apogee_rc_t *line)
{
        const double sin_phi = get_sin_phi(line);
        const double cos_phi = get_cos_phi(line);

        return line->pm_ra * cos_phi + line->pm_dec * sin_phi;
}

double mu_b_from_pa_dec_pm(const apogee_rc_t *line)
{
        const double sin_phi = get_sin_phi(line);
        const double cos_phi = get_cos_phi(line);

        return line->pm_dec * cos_phi - line->pm_ra * sin_phi;
}
#endif

double l_from_radec(const apogee_rc_t *line)
{
        return rad_to_deg(atan2(get_sin_phi_II(line), get_cos_phi_II(line))) + L_GAL;
}

double b_from_radec(const apogee_rc_t *line)
{
        return rad_to_deg(asin(get_sin_b(line)));
}

double mu_l_from_pa_dec_pm_II(const apogee_rc_t *line)
{
        double new_ra = line->ra + mas_to_rad(line->pm_ra);
        double new_dec = line->dec + mas_to_rad(line->pm_dec);
        apogee_rc_t new_line = *line;
        new_line.ra = new_ra;
        new_line.dec = new_dec;
        double res = deg_to_mas(l_from_radec(&new_line) - l_from_radec(line)) * cos(line->b);
        return res;
}

double mu_b_from_pa_dec_pm_II(const apogee_rc_t *line)
{
        double new_ra = line->ra + mas_to_rad(line->pm_ra);
        double new_dec = line->dec + mas_to_rad(line->pm_dec);
        apogee_rc_t new_line = *line;
        new_line.ra = new_ra;
        new_line.dec = new_dec;
        double res = deg_to_mas(b_from_radec(&new_line) - b_from_radec(line));
        return res;
}


static double delta_mu_b_by_ra_step(const apogee_rc_t *line)
{
	apogee_rc_t tmp = *line;
	tmp.ra += DRVT_STEP;

	return mu_b_from_pa_dec_pm_II(&tmp) - mu_b_from_pa_dec_pm_II(line);
}

static double delta_mu_b_by_dec_step(const apogee_rc_t *line)
{
	apogee_rc_t tmp = *line;
	tmp.dec += DRVT_STEP;

	return mu_b_from_pa_dec_pm_II(&tmp) - mu_b_from_pa_dec_pm_II(line);
}

static double delta_mu_l_by_ra_step(const apogee_rc_t *line)
{
	apogee_rc_t tmp = *line;
	tmp.ra += DRVT_STEP;

	return mu_l_from_pa_dec_pm_II(&tmp) - mu_l_from_pa_dec_pm_II(line);
}

static double delta_mu_l_by_dec_step(const apogee_rc_t *line)
{
	apogee_rc_t tmp = *line;
	tmp.dec += DRVT_STEP;

	return mu_l_from_pa_dec_pm_II(&tmp) - mu_l_from_pa_dec_pm_II(line);
}

/*
 * Error estimated translation
 */
double errors_ecliptic_to_gal_mu_l(const apogee_rc_t *line)
{
	const double drv_mu_l_ra = delta_mu_l_by_ra_step(line) / DRVT_STEP;
	const double drv_mu_l_dec = delta_mu_l_by_dec_step(line) / DRVT_STEP;
        return sqrt(pow_double(drv_mu_l_ra * line->pm_ra_err, 2) +
		    pow_double(drv_mu_l_dec * line->pm_dec_err, 2));
}

double errors_ecliptic_to_gal_mu_b(const apogee_rc_t *line)
{
	const double drv_mu_b_ra = delta_mu_b_by_ra_step(line) / DRVT_STEP;
	const double drv_mu_b_dec = delta_mu_b_by_dec_step(line) / DRVT_STEP;
        return sqrt(pow_double(drv_mu_b_ra * line->pm_ra_err, 2) +
		    pow_double(drv_mu_b_dec * line->pm_dec_err, 2));
}
