#include <stdlib.h>
#include <string.h>

#include "mem.h"
#define MAX_NAME_LENGH          16

char *generic_name(void)
{
        char *s = (char *)dv_alloc(sizeof(char) * MAX_NAME_LENGH);
        return s;
}
