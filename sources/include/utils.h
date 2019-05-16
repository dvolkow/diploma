#ifndef UTILS_H
#define UTILS_H         1

#include <stdbool.h>
#include "types.h"
#include "math.h"
#include "filter.h"

typedef struct {
        double theta;
        double err;
        double sd;
        unsigned int size;
} average_res_t;


void sort_iteration_storage_by_r(iteration_storage_t *,
                                 const size_t);

iteration_storage_t *iteration_storage_create(const apogee_rc_table_t *,
                                                const opt_t *);

apogee_rc_table_t *get_limited_generic(const void *,
                                       const filter_t *,
                                       filter_mode_t);

average_res_t *get_average_theta(const iteration_storage_t *,
                                 const opt_t *,
                                 const unsigned int);
filter_t *filter_factory(const parser_t *);

void update_table_R0(apogee_rc_table_t *, double);

#endif // UTILS_H
