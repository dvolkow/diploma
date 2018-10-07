#include <stdio.h>
#include <math.h>
#include <string.h>

#include "core.h"
#include "jtest.h"
#include "types.h"
#include "io.h"
#include "mem.h"
#include "math.h"
#include "debug.h"

double get_beta_n(const apogee_rc_t *line, beta_ord_t type)
{
        switch (type) {
                case FIRST:
                        return -cos(line->l) * cos(line->b);
                case SECOND:
                        return -sin(line->l) * cos(line->b);
                case THIRD:
                        return -sin(line->b);
                default:
                        printf("%s: type error!\n", __FUNCTION__);
        }

        return 0;
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

void fill_mnk_matrix_vr(double *matrix, const int ord,
                         apogee_rc_table_t *table)
{
        unsigned int i, j, k;
        unsigned int len = BETA_QTY + ord;
        matrix_line_t *line = dv_alloc(sizeof(matrix_line_t));

        for (j = 0; j < table->size; ++j) {
                for (i = 0; i < BETA_QTY; ++i) {
                        line->_[i] = get_beta_n(table->data + j, i);
                }

                for (i = BETA_QTY; i < ord; ++i) {
                        line->_[i] = get_alpha_n(table->data + j, table->r_0, i);
                }

                for (i = 0; i < len; ++i) {
                        double m = line->_[i];
                        for (k = 0; k < len; ++k) {
                                matrix[i * len + k] += line->_[k] * m;
                        }
                }
        }
#ifdef DEBUG
        print_matrix(matrix, len);
#endif
}

void fill_mnk_matrix(double *matrix, const int ord, 
                        apogee_rc_table_t *table, eq_mode_t mode)
{
        switch (mode) {
                case VR_MODE:
                        fill_mnk_matrix_vr(matrix, ord, table);
                        break;
                default:
                        printf("%s: mode not implemented!\n", __FUNCTION__);
                        exit(1);
        }
}



void run_all_jtests()
{
        UTEST_LINE_PRINT();
        printf("%s: start unit tests...\n", __FUNCTION__);

        unsigned int i;
        FORALL_JTEST_TABLE(jtest_table, i) {
                TRY_TEST(&jtest_table[i]);
        }
        printf("%s: unit test completed!\n", __FUNCTION__);
        UTEST_LINE_PRINT();
}

void get_solution(int argc, char *argv[])
{
        int size = atoi(argv[ORD_ARG]);
        double *matrix = dv_alloc(sizeof(double) * (size + BETA_QTY) * 
                                                   (size + BETA_QTY));
        apogee_rc_table_t *table = read_table(argv[INPUT_FILE_ARG]);
        if (table == NULL) {
                printf("%s: fail to open %s!\n",
                                __FUNCTION__, argv[INPUT_FILE_ARG]);
                return;
        }

#ifdef DEBUG
        printf("%s: ord = %d\n", __FUNCTION__, size);
#endif
        fill_mnk_matrix(matrix, size, table, VR_MODE);
}

int main(int argc, char *argv[])
{
        run_all_jtests();
        get_solution(argc, argv);
        return 0;
}
