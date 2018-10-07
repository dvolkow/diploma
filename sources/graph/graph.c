#include <stdio.h>
#include <math.h>

#include "graph.h"
#include "io.h"
#include "opt.h"
#include "math.h"


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
                PRINT_IO_OPEN_ERROR(OUTPUT_RESULT_FILENAME);
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


/*
static double get_sun_shift(const apogee_rc_table_t *table, 
                                const opt_t *solution, const int idx)
{
}
*/

void dump_rotation_curve(apogee_rc_table_t *table, opt_t *solution)
{
        FILE *fout = fopen(RC_OUT_FILE_NAME, "w");
        if (fout == NULL) {
                PRINT_IO_OPEN_ERROR(RC_OUT_FILE_NAME);
                return;
        }

        FILE *oout = fopen(DF_OUT_FILE_NAME, "w");
        if (oout == NULL) {
                PRINT_IO_OPEN_ERROR(RC_OUT_FILE_NAME);
                return;
        }

        FILE *sout = fopen(SUN_POINT_FILE_NAME, "w");
        if (sout == NULL) {
                PRINT_IO_OPEN_ERROR(SUN_POINT_FILE_NAME);
                return;
        }

        printf("%s: enrty\n", __func__);
        unsigned int i;
        double r;
        for (i = 0; i < table->size; ++i) {
                r = get_R_distance(&table->data[i], solution->r_0);
                double theta = r * ((table->data[i].v_helio - 
                                     solution->s.data[U_P] * get_beta_n(&table->data[i], FIRST) 
                                    + solution->s.data[W_P] * sin(table->data[i].b)) / 
                                        (solution->r_0 * sin(table->data[i].l) *
                                                             cos(table->data[i].b)) + OMEGA_SUN);
                fprintf(oout, "%lf %lf\n", r, theta);
        }

        r = ROTC_LOWER_BOUND;
        while (r < ROTC_UPPER_BOUND) {
                double theta = r * (OMEGA_SUN - solution->s.data[V_P] / solution->r_0) -
                                        2 * solution->s.data[A_P] * (r - solution->r_0);
                for (i = BETA_QTY + 1; i < solution->s.size; ++i) {
                        theta += solution->s.data[i] / dv_factorial(i - BETA_QTY + 1) * 
                                        pow_double(r - solution->r_0, i - BETA_QTY + 1);
                } 

                fprintf(fout, "%lf %lf\n", r, theta);
                r += ROTC_STEP_R;
        }

        fprintf(sout, "%lf %lf\n", solution->r_0, solution->r_0 * OMEGA_SUN);

        fclose(sout);
        fclose(oout);
        fclose(fout);
}
