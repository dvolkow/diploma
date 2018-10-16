#include "mem.h"
#include "asserts.h"
#include <stdlib.h>
#include <stdio.h>


static unsigned char __g_memory_buffer[MAX_MEMORY_USAGE];
static size_t current_usage = 0;



void *dv_alloc(const size_t size) 
{
        MEMORY_LIMIT_ASSERT(current_usage + size);
        MEMORY_ALLOC_ASSERT(size);
        void *ret = (void*)(__g_memory_buffer + current_usage);
        current_usage += size;
#ifdef DEBUG
        printf("%s: call from [%p], now usage %lu bytes\n",
                        __func__, __builtin_return_address(0),
                        current_usage);
#endif
        return ret;
}


#define G_START_SHIFT 8
void *dv_mm_get_current_top(void)
{
        return (void *)&__g_memory_buffer[current_usage > G_START_SHIFT 
                        ? current_usage - G_START_SHIFT 
                        : current_usage];
}

void dump_memory_usage(void)
{
        printf("%s: [debug]: Total usage %u byte, of %u (%0.3lf%)\n",
                        __func__, current_usage, MAX_MEMORY_USAGE,
                                current_usage / MAX_MEMORY_USAGE * 100);
}
