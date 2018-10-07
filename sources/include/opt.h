#ifndef OPT_H
#define OPT_H   1

#include "math.h"
#include "types.h"

typedef enum {
        LOWER,
        UPPER
} bound_t;

opt_t *opt_linear(linear_equation_t *eq, 
                                apogee_rc_table_t *table);

double lower_bound_search(linear_equation_t *eq, 
                                apogee_rc_table_t *table,
                                double r_0);
double upper_bound_search(linear_equation_t *eq, 
                                apogee_rc_table_t *table,
                                double r_0);
void get_errors(opt_t *, apogee_rc_table_t *table);
#endif // OPT_H
