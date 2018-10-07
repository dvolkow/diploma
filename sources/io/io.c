#include <stdio.h>
#include <stdlib.h>

#include "io.h"
#include "debug.h"


apogee_rc_table_t *read_table(const char *input_file_name)
{
        FILE *inp_f;
        unsigned int size = countlines(input_file_name);
        inp_f = fopen(input_file_name, "r");
        if (inp_f == NULL)
                return inp_f;

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
        }

        return table;
}

unsigned int countlines(char *filename) 
{
        FILE *f;
        f = fopen(filename, "a+");
        int ch = 0;
        unsigned int lines = 0;

        if (f == NULL)
                return 0;

        while ((ch = fgetc(f)) != EOF)
        {
                if (ch == '\n')
                        lines++;
        }
        fclose(f);
        return lines;
}

void output_result()
{
        return;
}
