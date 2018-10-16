#include <stdio.h>
#include <stdlib.h>

#include "io.h"
#include "math.h"
#include "mem.h"
#include "types.h"
#include "debug.h"


apogee_rc_table_t *read_table(const char *input_file_name)
{
        FILE *inp_f;
        unsigned int size = countlines(input_file_name);
        inp_f = fopen(input_file_name, "r");
        if (inp_f == NULL) {
                PRINT_IO_OPEN_ERROR(input_file_name);
                return NULL;
        }

        apogee_rc_t *apogee_rc = dv_alloc(sizeof(apogee_rc_t) * size);
        apogee_rc_table_t *table = dv_alloc(sizeof(apogee_rc_table_t));
        table->data = apogee_rc;
        table->size = size;

        unsigned int i;
        for (i = 0; i < size; ++i) {
                fscanf(inp_f, "%lf %lf %lf %lf %lf %lf",
                                &(apogee_rc[i].l),
                                &(apogee_rc[i].b),
                                &(apogee_rc[i].v_helio),
                                &(apogee_rc[i].dist),
                                &(apogee_rc[i].pm_ra),
                                &(apogee_rc[i].pm_dec)
                                );
                apogee_rc[i].l = deg_to_rad(apogee_rc[i].l);
                apogee_rc[i].b = deg_to_rad(apogee_rc[i].b);
                apogee_rc[i].id = i;
        }

#ifdef DEBUG
        //print_table(table);
        printf("%s: size table: %lu\n", __func__, table->size);
#endif
        return table;
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

        fscanf(fin, "%u", &solution->s.size); // 3 + ord!
        solution->s.data = dv_alloc(sizeof(double) * solution->s.size);

        fscanf(fin, "%lf", &solution->r_0); 
        fscanf(fin, "%lf", &solution->sq); 
        unsigned int i;
        for (i = 0; i < solution->s.size; ++i) {
                fscanf(fin, "%lf", &solution->s.data[i]); 
        }

        fscanf(fin, "%u", &solution->size); // size of sample
        return solution;
}

void output_result()
{
        return;
}
