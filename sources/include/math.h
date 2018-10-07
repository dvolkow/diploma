#ifndef MATH_H  
#define MATH_H  1

#include "types.h"

#define MAX_EQUATION_SIZE       10

typedef enum {
        EQUATION,
        SOLUTION,
} linear_type_t;

typedef struct {
        double  *data;
        int     size;
} linear_equation_t;

typedef struct {
        double  *data;        
        int     size;
} linear_eq_solve_t;

void solve(linear_equation_t *, linear_eq_solve_t *, linear_eq_solve_t *);

double get_R_distance(apogee_rc_t *line, double r_0);

#define PRECACHED_FACTORIAL_LEN         12
int dv_factorial(const int n);

double dot_prod(double *a, double *b, int size);

#endif // MATH_H

