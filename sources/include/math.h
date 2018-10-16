#ifndef MATH_H  
#define MATH_H  1

#include "types.h"
#include <stdlib.h>

#define M_PI 3.14159265358979323846


static inline double deg_to_rad(const double deg)
{
        return deg * M_PI / 180;
}

static inline double rad_to_deg(const double rad)
{
        return rad * 180 / M_PI;
}

static inline double pow_double(const double m, const unsigned int n) 
{
        unsigned int i;
        double _m = m;
        for (i = 0; i < n - 1; ++i) 
                _m *= m;
        return _m;
}



typedef enum {
        EQUATION,
        SOLUTION,
} linear_type_t;


typedef struct {
        double  *data;
        double  *right;
        int     size;
        short   ord;
} linear_equation_t;


typedef struct {
        double  *data;        
        int     size;
} linear_eq_solve_t;


typedef struct {
        double l;
        double h;
} prec_t;


typedef struct {
        linear_eq_solve_t s;
        double r_0;
        double sq;
        prec_t *bounds;
        size_t size;
} opt_t;


typedef struct {
        double x;
        double y;
        double z;
} point_t;


void solve(linear_equation_t *, linear_eq_solve_t *);
void inverse_and_diag(linear_equation_t *eq, linear_equation_t *res);

double get_R_distance(const apogee_rc_t *line, double r_0);
double get_error_mnk_estimated(const double p, __attribute__((__unused__)) const int nfree,
                                const double sd);
point_t *get_point(const apogee_rc_t *line); 

double get_median(const double *data, const size_t size);

double get_limit_by_eps(const unsigned int size);

#define PRECACHED_FACTORIAL_LEN         25
double dv_factorial(const int n);

double dot_prod(double *a, double *b, int size);

int math_init(void);
void math_exit(void);

#endif // MATH_H

