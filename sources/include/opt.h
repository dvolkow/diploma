#ifndef OPT_H
#define OPT_H   1

#include "math.h"
#include "types.h"

opt_t *opt_linear(linear_equation_t *eq, 
                                apogee_rc_table_t *table);

#endif // OPT_H
