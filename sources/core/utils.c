#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "types.h"
#include "mem.h"
#include "graph.h"
#include "math.h"
#include "utils.h"

double get_beta_n(const apogee_rc_t *line, beta_ord_t type);

int r_comparator(const void *p1, const void *p2)
{
        const iteration_storage_t *i1 = p1;
        const iteration_storage_t *i2 = p2;

        if (i1->r < i2->r)
                return -1;
        else if (i1->r > i2->r)
                return 1;
        else 
                return 0;
}

void sort_iteration_storage_by_r(iteration_storage_t *storage, const size_t size) 
{
        qsort(storage, size, sizeof(iteration_storage_t), r_comparator);
}


iteration_storage_t *iteration_storage_create(const apogee_rc_table_t *table, 
                                                const opt_t *solution)
{
        iteration_storage_t *storage = dv_alloc(sizeof(iteration_storage_t) *
                                        table->size);
        unsigned int i;
        for (i = 0; i < table->size; ++i) {
                double r = get_R_distance(&table->data[i], solution->r_0);
                storage[i].data = table->data[i];
                storage[i].r = r;
                storage[i].theta = r * ((table->data[i].v_helio - 
                                     solution->s.data[U_P] * get_beta_n(&table->data[i], FIRST) 
                                    + solution->s.data[W_P] * sin(table->data[i].b)) / 
                                        (solution->r_0 * sin(table->data[i].l) *
                                                             cos(table->data[i].b)) + OMEGA_SUN);
        }

        return storage;
}

average_res_t *get_average_theta(const iteration_storage_t *st_part, 
                                 const opt_t *solution,
                                 const size_t size) 
{
        assert(size != 0);

        double c = 0;
        double a = 0;
        unsigned int i;
        for (i = 0; i < size; ++i) {
                c += st_part[i].data.v_helio + 
                        solution->s.data[U_P] * cos(st_part[i].data.l) * cos(st_part[i].data.b)
                        + solution->s.data[W_P] * sin(st_part[i].data.b) +
                                        OMEGA_SUN * solution->r_0 * sin(st_part[i].data.l) *
                                                        cos(st_part[i].data.b);
                a += sin(st_part[i].data.l) * cos(st_part[i].data.b) * solution->r_0 / st_part[i].r;
        }
        average_res_t *res = dv_alloc(sizeof(average_res_t));
        res->theta = c / a;
        // TODO: errors fill
        res->err = 0;
        res->size = size;

        return res;
}
