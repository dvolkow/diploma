#include <stdio.h>
#include <math.h>

#include "core.h"
#include "jtest.h"
#include "types.h"
#include "io.h"
#include "mem.h"
#include "math.h"

#define UTEST_LINE              \
        "--------------------------------"

#define UTEST_LINE_PRINT()        \
        printf("%s\n", UTEST_LINE)


beta_coeff_t *get_beta_for_line(const apogee_rc_t *line)
{
        beta_coeff_t *beta = dv_alloc(sizeof(beta_coeff_t));

        beta->_1 = -cos(line->l) * cos(line->b);
        beta->_2 = -sin(line->l) * cos(line->b);
        beta->_3 = -sin(line->b);

        return beta;
}

double get_alpha_n(const apogee_rc_t *line, 
                   const double r_0,
                   const int n)
{
        double R = get_R_distance(line, r_0);
        int i;
        double res = sin(line->l) * cos(line->b) * r_0 / R * (R - r_0);
        if (n == 0)
                return res * (-2);

        for (i = 1; i < n; ++i) 
                res *= (R - r_0) / (i + 1);

        return res;
}

void run_all_jtests()
{
        UTEST_LINE_PRINT();
        printf("%s: start unit tests...\n", __FUNCTION__);
        jmemory();
        jmath();
        jio();
        apogee_rc_t *table = read_table(INPUT_TABLE_FILE_NAME);
        jtrigonometry(table);
        printf("%s: unit test completed!\n", __FUNCTION__);
        UTEST_LINE_PRINT();
}

int main(int argc, char *argv[])
{
        run_all_jtests();
        return 0;
}
