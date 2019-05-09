#ifndef MEM_H
#define MEM_H   1
#include <stdlib.h>
#include <string.h>

#define MAX_MEMORY_USAGE        150000000

void *dv_alloc(const size_t size);
void *dv_dalloc(const size_t size, const size_t num);
void *dv_mm_get_current_top(void);
void dump_memory_usage(void);

#endif // MEM_H
