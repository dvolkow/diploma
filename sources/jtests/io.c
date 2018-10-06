#include <stdio.h>
#include <stdlib.h>

#include "io.h"
#include "types.h"
#include "debug.h"

#define PRINT_JOI_RES(s)                        \
        printf("%s: count lines in %s: %d\n",   \
                        __FUNCTION__, (s), countlines((s)))

int jio()
{
        PRINT_JOI_RES("../sources/jtests/memory.c");
        PRINT_JOI_RES("./deleted_file");
        PRINT_JOI_RES("../sources/jtests/io_test.txt");

        apogee_rc_t *table = read_table(INPUT_TABLE_FILE_NAME);
        print_table(table, countlines(INPUT_TABLE_FILE_NAME));
}
