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
        assert(dist > 0);
        assert(r_0 > 0);
        double root = dist * dist * cos(line->b) * cos(line->b) +
                r_0 * r_0 - 2 * r_0 * dist * cos(line->l) * cos(line->b);

        assert(root >= 0);
        return sqrt(root);
}

static double get_x(const apogee_rc_t *line) 
{
        return line->dist * cos(line->l) * cos(line->b);
}

static double get_y(const apogee_rc_t *line) 
{
        return line->dist * sin(line->l) * cos(line->b);
}

static double get_z(const apogee_rc_t *line) 
{
        return line->dist * sin(line->b);
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

/*
 * Conversion between galactical and ecliptical velocities:
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

/*
 * Error estimated translation
 */
double errors_ecliptic_to_gal(const apogee_rc_t *line)
{
        const double ra_imp = pow_double(cos(deg_to_rad(DELTA_N)) * cos(line->ra - deg_to_rad(ALPHA_N)) / cos(line->b), 2);
        const double dec_imp = pow_double((sin(deg_to_rad(DELTA_N)) * sin(line->dec) / (cos(line->b) * cos(line->b)) - tan(line->b)) / (cos(line->dec) * cos(line->dec)), 2);

        return sqrt(ra_imp * pow_double(line->pm_ra_err, 2) + dec_imp * pow_double(line->pm_dec_err, 2));
}
