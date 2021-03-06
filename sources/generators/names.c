#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "mem.h"
#define MAX_NAME_LENGH          32

static char *generic_name(void)
{
        char *s = (char *)dv_alloc(sizeof(char) * MAX_NAME_LENGH);
        return s;
}

char *name_for_obj(const unsigned int i,
                   const unsigned int n,
                   const char *prefix)
{
        char *s = generic_name();
        sprintf(s, "%s_%u_%u", prefix, i, n);
        return s;
}

const char *basenames[] = {
        "N",
        "R_0",
        "SD",
        "U",
        "V",
        "W",
        "omega_0",
        "A"
};

char *name_for_column_uni_solution(const unsigned int i)
{
        char *s = generic_name();
        sprintf(s, "Theta%u", i);
        return s;
}
