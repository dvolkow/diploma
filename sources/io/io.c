#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "io.h"
#include "math.h"
#include "trigonometry.h"
#include "mem.h"
#include "types.h"
#include "debug.h"
#include "generators.h"

#include "core.h"
//#include "core_b.h"
//#include "core_l.h"
//#include "core_vr.h"


/**
 * TODO: update this
 * Table format:
 *
 * RA - DEC - GLON - GLAT - NVIS - VHELIO - RC_DIST - PM_RA - PM_DEC - PMMATCH
 */
apogee_rc_table_t *read_table(const char *input_file_name)
{
        FILE *inp_f;
        unsigned int size = countlines(input_file_name);
        inp_f = fopen(input_file_name, "r");
        if (inp_f == NULL) {
                PRINT_IO_OPEN_ERROR(input_file_name);
                return NULL;
        }

        apogee_rc_table_t *table = create_apogee_rc_table_by_size(size);
        apogee_rc_t *apogee_rc = table->data;

        unsigned int i;
        for (i = 0; i < size; ++i) {
//		TODO: double ra, dec, pm_ra, pm_dec, pm_ra_err, pm_dec_err;
                fscanf(inp_f, "%lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %d",
                                &(apogee_rc[i].ra),
                                &(apogee_rc[i].dec),
                                &(apogee_rc[i].l),
                                &(apogee_rc[i].b),
                                &(apogee_rc[i].nvis),
                                &(apogee_rc[i].v_helio),
                                &(apogee_rc[i].v_err),
                                &(apogee_rc[i].dist),
                                &(apogee_rc[i].pm_ra),
                                &(apogee_rc[i].pm_dec),
                                &(apogee_rc[i].pm_ra_err),
                                &(apogee_rc[i].pm_dec_err),
                                &(apogee_rc[i].pm_match)
                                );
#ifdef CALIBRATION
                apogee_rc[i].dist *= K_A;
#endif
                apogee_rc[i].l = deg_to_rad(apogee_rc[i].l);
                apogee_rc[i].b = deg_to_rad(apogee_rc[i].b);
                apogee_rc[i].ra = deg_to_rad(apogee_rc[i].ra);
                apogee_rc[i].dec = deg_to_rad(apogee_rc[i].dec);
                apogee_rc[i].pm_ra = apogee_rc[i].pm_ra / cos(apogee_rc[i].dec);

                apogee_rc[i].sin_l = sin(apogee_rc[i].l);
                apogee_rc[i].cos_l = cos(apogee_rc[i].l);
                apogee_rc[i].cos_b = cos(apogee_rc[i].b);
                apogee_rc[i].sin_b = sin(apogee_rc[i].b);
                apogee_rc[i].pm_l = K_PM * mu_l_from_pa_dec_pm_II(apogee_rc + i) * apogee_rc[i].cos_b;
                apogee_rc[i].pm_b = K_PM * mu_b_from_pa_dec_pm_II(apogee_rc + i);
                apogee_rc[i].pm_l_err = errors_ecliptic_to_gal_mu_l(apogee_rc + i);
                apogee_rc[i].pm_b_err = errors_ecliptic_to_gal_mu_b(apogee_rc + i);
#ifdef DEBUG
                printf("%u (%d): pm_l_I = %lf pm_b_I = %lf | pm_l_II = %lf pm_b_II = %lf\n", i,
                                apogee_rc[i].pm_match,
                                mu_l_from_pa_dec_pm(apogee_rc + i),
                                mu_b_from_pa_dec_pm(apogee_rc + i),
                                mu_l_from_pa_dec_pm_II(apogee_rc + i),
                                mu_b_from_pa_dec_pm_II(apogee_rc + i)
                                );
                printf("%u: l_cat = %lf b_cat = %lf | l = %lf b = %lf\n", i,
                                rad_to_deg(apogee_rc[i].l),
                                rad_to_deg(apogee_rc[i].b),
                                l_from_radec(&apogee_rc[i]),
                                b_from_radec(&apogee_rc[i])
                      );
                printf("%u: a_err = %lf d_err = %lf | l_err = %lf b_err = %lf | SD %lf vs %lf\n", i,
                                apogee_rc[i].pm_ra_err,
                                apogee_rc[i].pm_dec_err,
                                apogee_rc[i].pm_l_err,
                                apogee_rc[i].pm_b_err,
				sqrt(pow_double(apogee_rc[i].pm_ra_err, 2) + pow_double(apogee_rc[i].pm_dec_err, 2)),
				sqrt(pow_double(apogee_rc[i].pm_l_err, 2) + pow_double(apogee_rc[i].pm_b_err, 2))
                      );
#endif
        }

        return table;
}


void dump_table(const apogee_rc_table_t *table)
{
        parser_t *cfg = get_parser();
        const char *fname = cfg->dump_file_name;
        if (!strcmp(fname, cfg->input_file_name)) {
		PR_ERR("input and dump files are equals!");
		return;
        }

	if (fname == NULL) {
		PR_ERR("dump file option is empty! Use -d (--dump) arg.");
		return;
	}


        FILE *fout = fopen(fname, "w");
        CHECK_FILE_AND_RET(fout, fname);

        dsize_t i;
        for (i = 0; i < table->size; ++i) {
                fprintf(fout, "%.10e %.10e %.10e %.10e %d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %d\n",
                                rad_to_deg(table->data[i].ra),
                                rad_to_deg(table->data[i].dec),
                                rad_to_deg(table->data[i].l),
                                rad_to_deg(table->data[i].b),
                                table->data[i].nvis,
                                table->data[i].v_helio,
                                table->data[i].v_err,
                                table->data[i].dist,
				table->data[i].pm_ra * cos(table->data[i].dec),
				table->data[i].pm_dec,
				table->data[i].pm_ra_err,
				table->data[i].pm_dec_err,
				table->data[i].pm_match);
        }

        fclose(fout);
}



unsigned int countlines(const char *filename)
{
        FILE *f;
        f = fopen(filename, "a+");
        int ch = 0;
        unsigned int lines = 0;

        if (f == NULL) {
                PRINT_IO_OPEN_ERROR(filename);
                return 0;
        }

        while ((ch = fgetc(f)) != EOF)
        {
                if (ch == '\n')
                        lines++;
        }
        fclose(f);
        return lines;
}

/**
 * Solution file format:
 *      uint    s.size
 *      double  R_0
 *      double  sd(V_r)
 *      double . . .   (size params)
 *      . . .
 *      uint    size
 * See at sample.txt
 */
opt_t *read_solution(const char *input_file_name)
{
        FILE *fin = fopen(input_file_name, "r");
        if (fin == NULL) {
                PRINT_IO_OPEN_ERROR(input_file_name);
                return NULL;
        }

        opt_t *solution = dv_alloc(sizeof(opt_t));

        fscanf(fin, "%d", &solution->s.size); // 3 + ord!
        solution->s.data = dv_alloc(sizeof(double) * solution->s.size);

        fscanf(fin, "%lf", &solution->r_0);
        fscanf(fin, "%lf", &solution->sq);
        unsigned int i;
        for (i = 0; i < solution->s.size; ++i) {
                fscanf(fin, "%lf", &solution->s.data[i]);
        }

        fscanf(fin, "%u", &solution->size);
        return solution;
}
