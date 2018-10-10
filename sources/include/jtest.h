#ifndef JTESTS_H
#define JTESTS_H        1

#include "types.h"

#define UTEST_LINE              \
        "--------------------------------"

#define UTEST_LINE_PRINT()        \
        printf("%s\n", UTEST_LINE)


typedef struct {
        int (*test)(void);
        char *name;
} jtest_table_t;

// 
int jmemory(void);
int jmath(void);
int jio(void);
int jtrigonometry(void);
int jcore(void);


#define FORALL_JTEST_TABLE(table, i)                            \
        for (i = 0; table[i].test != NULL; ++i)

#define TRY_TEST(table_line)                                    \
        do {                                                    \
                if ((table_line)->test())                       \
                        printf("%s: fail for [%s]\n",           \
                                __func__,                       \
                                (table_line)->name);            \
                else                                            \
                        printf("%s: success for [%s]\n",        \
                                        __func__,               \
                                        (table_line)->name);    \
        } while (0)

int jtest_init();
void jtest_exit();
#endif // JTESTS_H
