#ifndef TYPES_H  
#define TYPES_H  1

#include <stdlib.h>
#include <stddef.h>

#include "core.h"

#define MAX_ORDER_SOLUTION      42

typedef unsigned int dsize_t;

/**
 * TODO: description for fields
 * @l:          longtitude, give by degress and transform 
 *              to radians for using into dipmloma project
 *
 * @b:          latitude, as l also read and make radins 
 *
 * @v_helio:    average velocity (kmps) 
 * @dist:       distance (pc)
 * @pm_ra:      proper motion to right ascention
 * @pm_dec:     proper motion to declination
 *
 * @eps:        error for solution parameter for this line
 */
typedef struct {
        // readable data:
        double  l; 
#define __mem_1 offsetof(apogee_rc_t, l)
        double  b;
#define __mem_2 offsetof(apogee_rc_t, b)
        double  v_helio;
#define __mem_3 offsetof(apogee_rc_t, v_helio)
        double  dist;
#define __mem_4 offsetof(apogee_rc_t, dist)
        double  pm_ra;
#define __mem_5 offsetof(apogee_rc_t, pm_ra)
        double  pm_dec;
#define __mem_6 offsetof(apogee_rc_t, pm_dec)
        double  ra;
        double  dec;
        double  v_err;
        double  pm_ra_err;
        double  pm_dec_err;

        // precalc data:
        double  pm_l;
        double  pm_b;
        double  pm_l_err;
        double  pm_b_err;
        
        double  sin_l;
        double  cos_l;
        double  sin_b;
        double  cos_b;
        
        // solution data:
        double  eps;
        double  vsd[3];
#define __mem_7 offsetof(apogee_rc_t, eps)
        // dsize_t id;
        int     nvis;
        int     pm_match;
} apogee_rc_t;


#define get_param(n, addr)         \
        (*(double *)(addr + __mem_##n))


/**
 * TODO: Special container structure needed for it
 * TEMPORARY SOLUTION?
 */
typedef struct {
        apogee_rc_t *data;
        double r_0;
        double w_sun;
        double omega_0;
        double n_err; // natural errors ^ 2
        double sigma[3];
        unsigned int size;
} apogee_rc_table_t;


typedef struct {
        // may be idx ? that about performance?
        apogee_rc_t     data;

        double  err;
        double  r; // R
        double  theta;
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
