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

void sort_iteration_storage_by_r(iteration_storage_t *storage, const size_t size);
iteration_storage_t *iteration_storage_create(const apogee_rc_table_t *table, 
                                                const opt_t *solution);

apogee_rc_table_t *get_limited_generic(const void *table, 
                                       const filter_t *filter,
                                       filter_mode_t mode);

filter_t *filter_factory(const parser_t *cfg);
bool __limited_by_l(const void *line, const double l, const double h);
#endif // UTILS_H
