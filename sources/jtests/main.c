#include "jtest.h"

static jtest_table_t jtest_table[] = {
        { .test = jmemory,      .name = "memory test" }
      , { .test = jmath,        .name = "math test" }
      , { .test = jio,          .name = "input/output test" }
      , { .test = jtrigonometry, .name = "trigonometry test" }
      , { .test = jcore,        .name = "core test" }
      // MUST BE LAST:
      , { .test = NULL,         .name = NULL }
      
};

int jtest_init()
{
        UTEST_LINE_PRINT();
        printf("%s: start unit tests...\n", __func__);

        unsigned int i;
        FORALL_JTEST_TABLE(jtest_table, i) {
                TRY_TEST(&jtest_table[i]);
        }
        printf("%s: unit test completed!\n", __func__);
        UTEST_LINE_PRINT();
        return 0;
}

void jtest_exit() {}
