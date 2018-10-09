#ifndef TYPES_H  
#define TYPES_H  1

#include <stdlib.h>
#include "core.h"

#define MAX_ORDER_SOLUTION      42

/**
 * TODO: description for fields
 *
 */
typedef struct {
        double l;
        double b;
        double v_helio;
        double dist;
        double pm_ra;
        double pm_dec;
} apogee_rc_t;

typedef struct {
        apogee_rc_t *data;
        double r_0; 
        size_t size;
} apogee_rc_table_t;


typedef struct {
        double  err; 
        double  r; // R 
        double  theta;

        apogee_rc_t     data;
        // may be idx ? that about performance?
} iteration_storage_t;

typedef struct {
        double _[BETA_QTY + MAX_ORDER_SOLUTION];
} matrix_line_t;

typedef struct {
        char *name;
} parameter_t;


enum {
        U_P = 0,
        V_P,
        W_P,
        A_P,
};

#endif // TYPES_H
