#ifndef UTILS_H
#define UTILS_H         1

#include <stdbool.h>
#include "types.h"
#include "math.h"

typedef struct {
        double theta;
        double err;
        unsigned int size;
} average_res_t;

typedef enum {
        L_FILTER
      , B_FILTER
      , ERR_FILTER
} filter_mode_t;


/**
 * @l: lower bound for parameter
 * @h: upper bound for parameter
 * @f: oracle for line of table that 
 *     const void* parameter
 */
typedef struct {
        double l;
        double h;
        bool (*f)(const void *, const double, const double);
} filter_t;

void sort_iteration_storage_by_r(iteration_storage_t *storage, const size_t size);

apogee_rc_table_t *get_limited_generic(const void *table, 
                                       const filter_t *filter,
                                       filter_mode_t mode);

bool __limited_by_l(const void *line, const double l, const double h);
#endif // UTILS_H
