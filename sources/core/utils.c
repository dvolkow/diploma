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
#ifdef DEBUG
	printf("%s: size = %u\n", __func__, table->size);
#endif
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


/**
 * Average theta by R for graphic plotting. Use
 * LSE by each interval. More details into text.
 */
static inline double __get_c(const iteration_storage_t *st_part,
                                 const opt_t *solution)
{
        return st_part->data.v_helio +
                        solution->s.data[U_P] * cos(st_part->data.l) * cos(st_part->data.b)
                        + solution->s.data[W_P] * sin(st_part->data.b) +
                                        OMEGA_SUN * solution->r_0 * sin(st_part->data.l) *
                                                        cos(st_part->data.b);
}

static inline double __get_a(const iteration_storage_t *st_part,
                                const opt_t *solution)
{
        return sin(st_part->data.l) * cos(st_part->data.b) * solution->r_0 / st_part->r;
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
                c += __get_c(&st_part[i], solution);
                a += __get_a(&st_part[i], solution);
        }

        average_res_t *res = dv_alloc(sizeof(average_res_t));
        res->theta = c / a;
        res->sd = 0;

        for (i = 0; i < size; ++i) {
                double d = __get_c(&st_part[i], solution) - 
                                __get_a(&st_part[i], solution) * res->theta;
                res->sd += pow_double(d, 2);
        }

        res->sd = sqrt(res->sd / (size - 1));
        res->err = res->sd / sqrt(2. * (size - 1)); 

        res->size = size;

        return res;
}


/**
 * ------ Filters for dataset --------
 */

/**
 * Before call must be set __cfg->h (high global limit)
 * as sd by sample
 */
bool __limited_by_eps(const void *line, 
                      const double l, 
                      const double h)
{
        return get_param(7, line) / l < h;
}


bool __limited_by_l(const void *line, const double l, const double h)
{
        return get_param(1, line) < deg_to_rad(h) && get_param(1, line) >= deg_to_rad(l);
}

bool __limited_by_b(const void *line, const double l, const double h)
{
        return get_param(2, line) < deg_to_rad(h) && get_param(2, line) >= deg_to_rad(l);
}

filter_t *filter_factory(const parser_t *cfg)
{
        filter_t *filter = dv_alloc(sizeof(filter));
        filter->l = cfg->l;
        filter->h = cfg->h;

        switch (cfg->filter) {
                case L_FILTER:
                        filter->f = __limited_by_l;
                        break;
                case B_FILTER:
                        filter->f = __limited_by_b;
                        break;
                case ERR_FILTER:
                        filter->f = __limited_by_eps;
                        break;
                default:
                        filter->f = NULL;
        }

        return filter;
}

/**
 * Return transorm table!
 */
apogee_rc_table_t *get_limited_replace(const void *table, const filter_t *filter)
{
        apogee_rc_table_t *src = table;
        unsigned int i;
        unsigned int count = 0;
        // no filter check:
        if (filter->f == NULL)
                return src;

        for (i = 0; i < src->size; ++i) {
                if (filter->f(&src->data[i], filter->l, filter->h)) {
                        src->data[count++] = src->data[i];
                } 
        }

        src->size = count;
        return src;
}

apogee_rc_table_t *get_limited_generic(const void *table, 
                                       const filter_t *filter,
                                       filter_mode_t mode) 
{
        if (filter == NULL)
                return (apogee_rc_table_t *)table;

        switch (mode) {
                case L_FILTER:
                default:
                        return get_limited_replace(table, filter);
                        // TODO: mode supporting
                        printf("%s: unknown mode 0x%x!\n",
                                        __func__, mode);
                        return NULL;
        }
}
