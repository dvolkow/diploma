#include "math.h"
#include "types.h"

#include <math.h>

double get_R_distance(apogee_rc_t *line, double r_0)
{
        return sqrt(line->dist * line->dist * cos(line->b) *
                        cos(line->b) + r_0 * r_0 +
                        2 * r_0 * line->b * cos(line->l) * cos(line->b));
}


