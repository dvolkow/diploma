#ifndef MEM_H
#define MEM_H   1
#include <stdlib.h>

#define MAX_MEMORY_USAGE        150000000

void *dv_alloc(const size_t size);
void *dv_mm_get_current_top(void);

#endif // MEM_H
