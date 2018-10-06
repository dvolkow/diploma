#ifndef MATH_H  
#define MATH_H  1

#define MAX_EQUATION_SIZE       10

typedef enum {
        EQUATION,
        SOLUTION,
} linear_type_t;

typedef struct {
        double  data[MAX_EQUATION_SIZE];
        int     size;
} linear_equation_t;

typedef struct {
        double  data[MAX_EQUATION_SIZE];
        int     size;
} linear_eq_solve_t;

void solve(linear_equation_t *eq, linear_eq_solve_t *s);

#endif // MATH_H

