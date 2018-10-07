#ifndef UTILS_H
#define UTILS_H         1

#include "types.h"
#include "math.h"

typedef struct {
        double theta;
        double err;
        unsigned int size;
} average_res_t;

void sort_iteration_storage_by_r(iteration_storage_t *storage, const size_t size);
#endif // UTILS_H
