#include "math.h"
#include "mem.h"
#include <gsl/gsl_linalg.h>

void *make_linear_struct(double *data, int size, 
                                linear_type_t type)
{
        void *block;

        switch(type) {
        case EQUATION:
                block = dv_alloc(sizeof(linear_equation_t)); 
                break;
        case SOLUTION:
                block = dv_alloc(sizeof(linear_eq_solve_t)); 
                break;
        }
        return block;
}

void solve(linear_equation_t *eq, linear_eq_solve_t *s)
{

}
