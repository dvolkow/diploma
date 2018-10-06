#include <stdio.h>

#include "core.h"
#include "jtest.h"

#define UTEST_LINE              \
        "--------------------------------"

#define UTEST_LINE_PRINT()        \
        printf("%s\n", UTEST_LINE)


void run_all_jtests()
{
        UTEST_LINE_PRINT();
        printf("%s: start unit tests...\n", __FUNCTION__);
        jmemory();
        jmath();
        jio();
        printf("%s: unit test completed!\n", __FUNCTION__);
        UTEST_LINE_PRINT();
}

int main(int argc, char *argv[])
{
        run_all_jtests();
        return 0;
}
