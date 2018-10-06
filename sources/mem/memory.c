#include "mem.h"
#include "asserts.h"
#include <stdlib.h>


static unsigned char __g_memory_buffer[MAX_MEMORY_USAGE];
static size_t current_usage = 0;



void *dv_alloc(const size_t size) 
{
        MEMORY_LIMIT_ASSERT(current_usage + size);
        void *ret = (void*)(__g_memory_buffer + current_usage);
        current_usage += size;
        return ret;
}

