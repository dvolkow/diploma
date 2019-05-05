#include <assert.h>
#include "types.h"
#include "core.h"
#include "math.h"
#include "io.h"
#include "unicore.h"

double get_beta_n(const apogee_rc_t *line, beta_ord_t type);
double get_alpha_n(const apogee_rc_t *line, 
                   const double r_0,
                   const int n);

int jcore()
{
        apogee_rc_table_t *table = read_table(INPUT_TABLE_FILE_NAME);
        assert(table->size == 1);


        unsigned int i;
        for (i = 0; i < BETA_QTY; ++i)
                printf("%s: get_beta_n(%d) for line table default: %lf\n",
                        __func__, i, core_vr_get_beta_n(&table->data[0], i));
        for (i = 1; i < 10; ++i)
                printf("%s: get_alpha_n(%d) for line table default: %lf\n",
                        __func__, i, core_vr_get_alpha_n(&table->data[0], i, 8));
        return 0;
}
