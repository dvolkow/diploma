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
//        printf("old l = %lf, new l = %lf\n", l_from_radec(line), l_from_radec(&new_line));
//        printf("old l - new l = %lf\n", res);
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
//        printf("old b = %lf, new b = %lf\n", b_from_radec(line), b_from_radec(&new_line));
//        printf("old b - new b = %lf\n", b_from_radec(&new_line) - b_from_radec(line));
        return res;
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
