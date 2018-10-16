#include <math.h>
#include <assert.h>

#include "types.h"
#include "math.h"
#include "mem.h"

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
