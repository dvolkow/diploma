#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "core.h"
#include "jtest.h"
#include "types.h"
#include "io.h"
#include "mem.h"
#include "math.h"
#ifdef DEBUG
#include "debug.h"
#endif // DEBUG
#include "opt.h"

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

static double s_alpha_n(const double R, 
                        const double sinl, 
                        const double cosb, 
                        const double r_0, const int n)
{
        if (n == 1)
                return -2 * (R - r_0) * r_0 * sinl * cosb / R;
        else 
                return r_0 * pow_double(R - r_0, n) * sinl * cosb / (R * dv_factorial(n));
}

double get_alpha_n(const apogee_rc_t *line, 
                   const double r_0,
                   const int n)
{
        assert(n > 0);
        
        double R = get_R_distance(line, r_0);
        return s_alpha_n(R, sin(line->l), cos(line->b), r_0, n);
}

void fill_mnk_matrix_vr(linear_equation_t *eq,
                         apogee_rc_table_t *table)
{
        unsigned int i, j, k;
        unsigned int len = eq->size;
        matrix_line_t *line = (matrix_line_t *)dv_alloc(sizeof(double) * len);

        memset(eq->data, 0, sizeof(double) * eq->size * eq->size);
        memset(eq->right, 0, sizeof(double) * eq->size);

        for (j = 0; j < table->size; ++j) {
                for (i = 0; i < BETA_QTY; ++i) {
                        line->_[i] = get_beta_n(&table->data[j], i);
                }

                for (i = BETA_QTY; i < eq->size; ++i) {
                        line->_[i] = get_alpha_n(&table->data[j], table->r_0, i - BETA_QTY + 1);
                }

                for (i = 0; i < len; ++i) {
                        double m = line->_[i];
                        for (k = 0; k < len; ++k) {
                                eq->data[i * len + k] += line->_[k] * m;
                        }
                        eq->right[i] += table->data[j].v_helio * m;
                }
        }
#ifdef DEBUG
        print_matrix(eq->data, len);
#endif
}

void fill_mnk_matrix(linear_equation_t *eq, 
                        apogee_rc_table_t *table, eq_mode_t mode)
{
        switch (mode) {
                case VR_MODE:
                        fill_mnk_matrix_vr(eq, table);
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
        assert(table->size != 0);

        linear_equation_t eq = {
                .data = matrix,
                .right = dv_alloc(sizeof(double) * size + BETA_QTY),
                .size = size + BETA_QTY,
                .ord = size
        };
#ifdef DEBUG
        printf("%s: ord = %d\n", __FUNCTION__, eq.ord);
#endif
        opt_t *solution = opt_linear(&eq, table);
        printf("%s: solution: R_0 = %lf\n",
                        __FUNCTION__, solution->r_0);


        printf("%s: solution: -R_0 = %lf\n", 
                        __FUNCTION__, lower_bound_search(&eq, table, solution->r_0));
        printf("%s: solution: +R_0 = %lf\n", 
                        __FUNCTION__, upper_bound_search(&eq, table, solution->r_0));
}

int main(int argc, char *argv[])
{
        run_all_jtests();
        get_solution(argc, argv);
        return 0;
}
