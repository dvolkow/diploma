#include <math.h>
#include <assert.h>

#include "math.h"
#include "types.h"

#ifdef DEBUG
#include "debug.h"
#endif

double get_R_distance(apogee_rc_t *line, const double r_0)
{
        double dist = line->dist;
        assert(dist > 0);
        assert(r_0 > 0);
        double root = dist * dist * cos(line->b) * cos(line->b) +
                r_0 * r_0 - 2 * r_0 * dist * cos(line->l) * cos(line->b);
        assert(root >= 0);
        return sqrt(root);
}


