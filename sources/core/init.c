#include "init.h"
#include "math.h"
#include "jtest.h"
#include "generators.h"

static init_t g_init_deinit_table[] = {
        { "math subsystem",  math_init, math_exit }
      , { "random generators subsystem",  random_seed_init, random_seed_exit }
      , { "unit test subsystem",  jtest_init, jtest_exit }

        // MUST BE LAST:
      , { NULL, NULL, NULL }
};

static int g_init_deinit_table_size = 0;

static int initialization_process()
{
        int i;
        int res = 0;
        for (i = 0; g_init_deinit_table[i].name != NULL; ++i) {
                g_init_deinit_table_size++;
                res = g_init_deinit_table[i].init();
                if (res)
                        goto fail;
                printf("%s: initialization for %s success!\n",
                                __func__, g_init_deinit_table[i].name);
        }
        return 0;
fail:
        printf("%s: failure for [%s] (#%d)!\n",
                        __func__, g_init_deinit_table[i].name,
                        i);
        for (; i >= 0; --i)
                g_init_deinit_table[i].exit();

        return res;
}

static void deinitialization_process()
{
        int i = g_init_deinit_table_size - 1;
        for (; i >= 0; --i) {
                g_init_deinit_table[i].exit();
                printf("%s: deinitialization for %s success!\n",
                                __func__, g_init_deinit_table[i].name);
        }
}


int main(int argc, char *argv[])
{
        initialization_process();
#ifdef DEBUG
        run_all_jtests();
#endif
        //get_solution(argc, argv);
        deinitialization_process();
        return 0;
}
