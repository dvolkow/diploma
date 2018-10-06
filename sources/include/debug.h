#ifndef DEBUG_H  
#define DEBUG_H  1

#include <stdio.h>
#include "io.h"

#define PRINTOK()       \
        printf("%s: line %d OK!\n", __FUNCTION__, __LINE__)

static inline void print_table(const apogee_rc_t *table, 
                                const int size)
{
        int i = 0;
        if (table == NULL)
                return;

        for (i = 0; i < size; ++i)
                printf("%s: %lf, %lf, %lf, %lf, %lf, %lf\n",
                                __FUNCTION__,
                                table[i].l,
                                table[i].b,
                                table[i].v_helio,
                                table[i].dist,
                                table[i].pm_ra,
                                table[i].pm_dec
                                );
}


#endif // DEBUG_H
