#ifndef ASSERTS_H
#define ASSERTS_H       1

#include "mem.h"
#include "math.h"
#include <assert.h>


#define MEMORY_LIMIT_ASSERT(a)         \
        assert((a) < MAX_MEMORY_USAGE)
        
#define FACTORIAL_STORAGE_ASSERT(a)     \
        assert((a) < PRECACHED_FACTORIAL_LEN)

#define MEMORY_ALLOC_ASSERT(a)          \
        assert((a) != 0)
#endif // ASSERTS_H
