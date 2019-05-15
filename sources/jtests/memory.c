#include <stdio.h>
#include "mem.h"
#include "math.h"
#include "asserts.h"
#include "types.h"


int jmemory()
{
      linear_equation_t *jle = dv_alloc(sizeof(linear_equation_t));
      linear_eq_solve_t *jes = dv_alloc(sizeof(linear_eq_solve_t));
      linear_equation_t *jle2 = dv_alloc(sizeof(linear_equation_t));

      assert((void *)jle != (void *)jes);
      assert(jes != NULL);
      assert(jle != NULL);
      printf("%s: linear_equation_t has %lu\n",
                      __func__, sizeof(*jle));
      printf("%s: linear_eq_solve_t has %lu\n",
                      __func__, sizeof(*jes));
      printf("%s: jle has %p memaddr\n",
                      __func__, (void *)jle);
      printf("%s: jes has %p memaddr\n",
                      __func__, (void *)jes);
      printf("%s: jle2 has %p memaddr\n",
                      __func__, (void *)jle2);
      return 0;
}
