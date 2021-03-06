#include <assert.h>
#include "init.h"
#include "math.h"
#include "debug.h"
#include "jtest.h"
#include "generators.h"
#include "graph.h"
#include "core.h"
#include "core_l.h"
#include "core_b.h"
#include "io.h"
#include "unicore.h"
#include "db.h"

static init_t g_init_deinit_table[] = {
        { "math subsystem",  math_init, math_exit }
      , { "parser subsystem",  parser_init, parser_exit }
      , { "random generators subsystem",  random_seed_init, random_seed_exit }
      , { "data base",  db_init, db_exit }
#ifdef JTEST_DEBUG
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


int main(int argc, const char *argv[])
{
        initialization_process();
        parse_args(argc, (char **)argv);
        parser_t *cfg = get_parser();
        if (cfg == NULL || !parser_t_is_valid(cfg)) {
                return -1;
        }

        apogee_rc_table_t *table = read_table(cfg->input_file_name);
        if (table == NULL) {
                printf("%s: fail to open %s!\n",
                                __func__, cfg->input_file_name);
                return -2;
        }

        assert(table->size != 0);

        db_add(generic_table()); // ERROR_LIMITED

        solution_mode_t mode = GET_SOLUTION_MODE(cfg);

        switch(mode) {
	case VR_PART_MODE:
		get_partial_vr_solution(table);
		break;
        case L_PART_MODE:
		get_partial_l_solution(table);
		break;
        case B_PART_MODE:
		get_partial_b_solution(table);
		break;
        case UNI_MODE:
                get_united_solution(table);
                break;
        case UNINAT_MODE:
                get_united_sigma_0_solution(table);
		break;
        case FIND_SIGMA0_MODE:
                find_united_sigma_0_solution(table);
		break;
        default:
		printf("Use -s [--solution] option to set solution mode.\nExit.\n");
		break;
        }

        const apogee_rc_table_t *dumped = db_get(ERROR_LIMITED);
        dump_objects_xyz(dumped, dumped->size, "err_xyz.txt");
        dump_objects_xyz(table, table->size, "final_xyz.txt");
        dump_table(table);
        deinitialization_process();
#ifdef DEBUG
        dump_memory_usage();
#endif
        return 0;
}
