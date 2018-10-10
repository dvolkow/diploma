#include "init.h"
#include "math.h"
#include "generators.h"

static init_t g_init_deinit_table[] = {
        { "math sybsystem",  math_init, math_exit }
      , { "random generators sybsystem",  random_seed_init, random_seed_exit }
        // MUST BE LAST:
      , { NULL, NULL, NULL }
};

static int g_init_deinit_table_size = 0;

int initialization_process()
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
        for (; i >= 0; --i)
                g_init_deinit_table[i].exit();

        return res;
}

void deinitialization_process()
{
        int i = g_init_deinit_table_size - 1;
        for (; i >= 0; --i) {
                g_init_deinit_table[i].exit();
                printf("%s: deinitialization for %s success!\n",
                                __func__, g_init_deinit_table[i].name);
        }
}
