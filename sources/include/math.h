#ifndef MATH_H  
#define MATH_H  1

#include "types.h"

#define MAX_EQUATION_SIZE       10

#define M_PI 3.14159265358979323846


static inline double deg_to_rad(const double deg)
{
        return deg * M_PI / 180;
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
        linear_eq_solve_t s;
        double r_0;
        double sq;
} opt_t;

void solve(linear_equation_t *, linear_eq_solve_t *);

double get_R_distance(apogee_rc_t *line, double r_0);

#define PRECACHED_FACTORIAL_LEN         12
int dv_factorial(const int n);

double dot_prod(double *a, double *b, int size);

#endif // MATH_H

