#ifndef ASSERTS_H
#define ASSERTS_H       1

#include "mem.h"
#include <assert.h>


#define MEMORY_LIMIT_ASSERT(a)         \
        assert((a) < MAX_MEMORY_USAGE)
        

#endif // ASSERTS_H
