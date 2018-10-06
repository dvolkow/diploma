#ifndef TYPES_H  
#define TYPES_H  1

#include "core.h"

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
        double err[MAX_ORDER_SOLUTION]; 
        double r; // R 
        double theta;
} iteration_storage_t;


typedef struct {
        double _1;
        double _2;
        double _3;
} beta_coeff_t;

typedef struct {
        double _[MAX_ORDER_SOLUTION];
} alpha_coeff_t;

typedef struct {
        beta_coeff_t *beta;
        alpha_coeff_t *alpha;
} matrix_line_mnk_t;
#endif // TYPES_H
