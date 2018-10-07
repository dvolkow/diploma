#include <stdio.h>
#include <math.h>

#include "graph.h"
#include "io.h"


#define OUTPUT_RESULT_FILENAME  \
        "result.txt"

parameter_t g_ptable[] = {
        { "u" }
      , { "v" }
      , { "w" }
      , { "A" }
      , { "th2" }
      , { "th3" }
      , { "th4" }
      , { "th5" }
      , { "th6" }
      , { "th7" }
      , { "th8" }
      , { "th9" }
      , { "th10" }
        // MUST BE LAST:
      , { NULL }
};

void dump_result(opt_t *opt, apogee_rc_table_t *table,
                  prec_t *p) 
{
        FILE *fout = fopen(OUTPUT_RESULT_FILENAME, "w"); 
        if (fout == NULL) {
                printf("%s: fail to open %s!\n", 
                                __func__, OUTPUT_RESULT_FILENAME);
                return;
        }

        fprintf(fout, "Result for APOGEE-RC dataset by %lu obj:\n",
                        table->size);
        fprintf(fout, "Optimal R_0: %lf (+%lf -%lf)\n",
                        opt->r_0, p->h - opt->r_0, opt->r_0 - p->l);
        fprintf(fout, "SD V_R: %lf kmps\n",
                        sqrt(opt->sq / (table->size + opt->s.size + 1)));
        unsigned int i;
#ifdef DEBUG
        printf("%s: opt->s.size = %u\n",
                        __func__, opt->s.size);
        for (i = 0; i < opt->s.size; ++i) {
                printf("%s: [%u] = %lf\n", __func__, 
                                i, opt->s.data[i]);
        }
#endif
        for (i = 0; i < opt->s.size; ++i) {
                fprintf(fout, "%s: %f (pm %f) \n",
                        g_ptable[i].name,
                        opt->s.data[i],
                        opt->bounds[i].l);
        }
        fclose(fout);
}
