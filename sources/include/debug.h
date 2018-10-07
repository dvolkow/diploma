#ifndef DEBUG_H  
#define DEBUG_H  1

#include <stdio.h>
#include "types.h"

#define PRINTOK()       \
        printf("%s: line %d OK!\n", __FUNCTION__, __LINE__)

static inline void print_table(const apogee_rc_table_t *table) 
{
        int i = 0;
        if (table == NULL)
                return;

        for (i = 0; i < table->size; ++i)
                printf("%s: %lf, %lf, %lf, %lf, %lf, %lf\n",
                                __FUNCTION__,
                                table->data[i].l,
                                table->data[i].b,
                                table->data[i].v_helio,
                                table->data[i].dist,
                                table->data[i].pm_ra,
                                table->data[i].pm_dec
                                );
}

static inline print_matrix(const double *matrix, const int size)
{
        int i, j;
        for (i = 0; i < size; ++i) {
                printf("%s: line %d: ", __FUNCTION__, i);
                for (j = 0; j < size; ++j) {
                        printf(" %lf", matrix[i * size + j]);
                }
                printf("\n");
        }
}

#endif // DEBUG_H
