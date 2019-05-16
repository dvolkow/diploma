#ifndef MATH_H
#define MATH_H  1

#include "types.h"
#include <stdlib.h>
#include <math.h>

#define M_PI            3.14159265358979323846
#define YR_TO_SEC       3.154e7

#define DRVT_STEP	1e-6

typedef struct {
	double x;
	double y;
} xy_t;


static inline double deg_to_rad(const double deg)
{
        return deg * M_PI / 180;
}

static inline double mas_to_rad(const double mas)
{
        return mas * 0.001 * M_PI / (180 * 60 * 60);
}

static inline double mas_to_deg(const double mas)
{
        return mas / (60 * 60 * 1000);
}

static inline double deg_to_mas(const double deg)
{
        return deg * 60 * 60 * 1000;
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


#define DOUBLE_COMPARE_EPSILON          (1e-8)

#define DOUBLE_EQUAL_EPS(a, b)          \
        (fabs((a) - (b)) < DOUBLE_COMPARE_EPSILON)


typedef enum {
        EQUATION,
        SOLUTION,
} linear_type_t;


typedef struct {
        double  *data;
        double  *right;
        unsigned int size;
        unsigned int ord;
} linear_equation_t;


typedef struct {
        double  *data;
        unsigned int size;
} linear_eq_solve_t;


typedef struct {
        double l;
        double h;
} prec_t;


typedef struct {
        linear_eq_solve_t s;
        double r_0;
        double dr_0;
        double sq;
        prec_t *bounds;
        unsigned int size;
} opt_t;

#define GET_SOLUTION_R0(p_opt)  \
        ((p_opt)->r_0)

typedef struct {
        double x;
        double y;
        double z;
} point_t;


void solve(linear_equation_t *, linear_eq_solve_t *);
void inverse_and_diag(linear_equation_t *eq, linear_equation_t *res);

double get_R_distance(const apogee_rc_t *line, double r_0);
double get_error_mnk_estimated(const double p, __attribute__((__unused__)) const unsigned int nfree,
                                const double sd);
point_t *get_point(const apogee_rc_t *line);

double get_median(const double *data, const unsigned int size);
double get_sd(const double *data, const unsigned int size);
double get_mean(const double *data, const unsigned int size);

double get_limit_by_eps(const unsigned int size);

#define PRECACHED_FACTORIAL_LEN         25
double dv_factorial(const unsigned int n);

double dot_prod(double *a, double *b, int size);

void add_matrix_to_matrix(const linear_equation_t *src,
                          linear_equation_t *dst);
int math_init(void);
void math_exit(void);

void linear_interpolation(const xy_t *,
			  const xy_t *,
			  xy_t *);
#endif // MATH_H

