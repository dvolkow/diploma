#include "init.h"
#include "math.h"
#include "debug.h"
#include "jtest.h"
#include "generators.h"
#include "core.h"

static init_t g_init_deinit_table[] = {
        { "math subsystem",  math_init, math_exit }
      , { "parser subsystem",  parser_init, parser_exit }
      , { "random generators subsystem",  random_seed_init, random_seed_exit }
#ifdef DEBUG
      , { "unit test subsystem",  jtest_init, jtest_exit }
#endif

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
#ifdef DEBUG
                printf("%s: initialization for %s success!\n",
                                __func__, g_init_deinit_table[i].name);
#endif
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
#ifdef DEBUG
                printf("%s: deinitialization for %s success!\n",
                                __func__, g_init_deinit_table[i].name);
#endif
        }
}


int main(int argc, char *argv[])
{
        initialization_process();
        parse_args(argc, argv);
        parser_t *cfg = get_parser();
        if (cfg == NULL || !parser_t_is_valid(cfg)) {
                printf("%s: invalid parameters!\n",
                                __func__);
                return -1;
        }

        g_mode_t mode = GET_MODE(cfg);
        switch(mode) {
                case SIMPLE_MODE:
                case ITERATE_MODE:
                        get_solution();
                        break;
                case GENERATION_MODE:
                        generate();
                default:
                        break;

        }

        deinitialization_process();
#ifdef DEBUG
        dump_memory_usage();
#endif
        return 0;
}
