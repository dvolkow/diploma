#include "mem.h"
#include "math.h"
#include "asserts.h"
#include "types.h"


int jmemory()
{
      linear_equation_t *jle = dv_alloc(sizeof(linear_equation_t));  
      linear_eq_solve_t *jes = dv_alloc(sizeof(linear_eq_solve_t));
      linear_equation_t *jle2 = dv_alloc(sizeof(linear_equation_t));
      
      assert(jle != jes);
      assert(jes != NULL);
      assert(jle != NULL);
      printf("%s: linear_equation_t has %lu\n", 
                      __FUNCTION__, sizeof(*jle));
      printf("%s: linear_eq_solve_t has %lu\n", 
                      __FUNCTION__, sizeof(*jes));
      printf("%s: jle has %p memaddr\n",
                      __FUNCTION__, jle);
      printf("%s: jes has %p memaddr\n",
                      __FUNCTION__, jes);
      printf("%s: jle2 has %p memaddr\n",
                      __FUNCTION__, jle2);
      return 0;
}
