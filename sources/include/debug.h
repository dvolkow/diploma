#ifndef DEBUG_H
#define DEBUG_H  1

#include <stdio.h>
#include "types.h"
#include "math.h"
#define PRINTOK()					\
        printf("%s: line %d OK!\n", __func__, __LINE__)

#define PR_WARN(message)				\
	printf("%s[%d]: WARNING: %s\n", __func__,	\
					__LINE__,	\
					(message))

#define PR_ERR(message)					\
	printf("%s[%d]: ERROR: %s\n",   __func__,	\
					__LINE__,	\
					(message))


static inline void print_table(const apogee_rc_table_t *table)
{
        int i = 0;
        if (table == NULL)
                return;

        for (i = 0; i < table->size; ++i)
                printf("%s: %lf, %lf, %lf, %lf, %lf, %lf\n",
                                __func__,
                                table->data[i].l,
                                table->data[i].b,
                                table->data[i].v_helio,
                                table->data[i].dist,
                                table->data[i].pm_ra,
                                table->data[i].pm_dec
                                );
}

static inline void print_vector(const double *data,
                                const unsigned int size)
{
#if 0
        unsigned int i;
        for (i = 0; i < size; ++i) {
                printf("%s: [debug] %lf\n", __func__,
                                        data[i]);
        }
#endif
}

static inline void print_matrix(const double *matrix, const int size)
{
#if 0
        int i, j;
        for (i = 0; i < size; ++i) {
                printf("%s: line %d: ", __func__, i);
                for (j = 0; j < size; ++j) {
                        printf(" %lf", matrix[i * size + j]);
                }
                printf("\n");
        }
#endif
}


void dump_memory_usage(void);

#endif // DEBUG_H
