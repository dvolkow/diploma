#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#include "core.h"
#include "types.h"
#include "io.h"
#include "mem.h"
#include "math.h"
#include "debug.h"
#include "opt.h"
#include "graph.h"
#include "utils.h"
#include "generators.h"

static matrix_line_t g_matrix_line;

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
                        printf("%s: type error!\n", __func__);
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
#ifdef DEBUG
        assert(n > 0);
#endif
        
        double R = get_R_distance(line, r_0);
        return s_alpha_n(R, sin(line->l), cos(line->b), r_0, n);
}

void fill_mnk_matrix_vr(linear_equation_t *eq,
                         apogee_rc_table_t *table)
{
        unsigned int i, j, k;
        unsigned int len = eq->size;
        matrix_line_t *line = &g_matrix_line;
#ifdef DEBUG
        assert(len <= BETA_QTY + MAX_ORDER_SOLUTION);
#endif

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
                        printf("%s: mode not implemented!\n", __func__);
                        exit(1);
        }
}

static void filter_get_and_apply(apogee_rc_table_t *table)
{
        /* TODO: Filter assigner: */
        parser_t *cfg = get_parser();

        switch (cfg->filter) {
                case L_FILTER:
                case B_FILTER:
                case ERR_FILTER:
                        table = get_limited_generic(table, 
                                                    filter_factory(cfg), 
                                                    L_FILTER);
                        break;
                default:
                        return;
        }
}

static opt_t *__get_solution_iterate(apogee_rc_table_t *table, 
                                        linear_equation_t *eq)
{
        assert(table->size != 0);

        filter_get_and_apply(table);
        opt_t *solution = opt_linear(eq, table);

        prec_t p = {
                .l = lower_bound_search(eq, table, solution->r_0),
                .h = upper_bound_search(eq, table, solution->r_0)
        };

        get_errors(solution, table);
        iteration_storage_t *st = iteration_storage_create(table, solution);
        dump_all(solution, &p, st);
        return solution;
}


void get_solution()
{
        parser_t *cfg = get_parser();
        int size = cfg->ord;
        double *matrix = dv_alloc(sizeof(double) * (size + BETA_QTY) * 
                                                   (size + BETA_QTY));
        apogee_rc_table_t *table = read_table(cfg->input_file_name);


        if (table == NULL) {
                printf("%s: fail to open %s!\n",
                                __func__, cfg->input_file_name);
                return;
        }
        assert(table->size != 0);
        
        filter_get_and_apply(table);

        dump_objects_xyz(table, table->size);
        dump_table(table);

        linear_equation_t eq = {
                .data = matrix,
                .right = dv_alloc(sizeof(double) * size + BETA_QTY),
                .size = size + BETA_QTY,
                .ord = size
        };
        opt_t *solution = opt_linear(&eq, table);
#ifdef DEBUG
        printf("%s: solution: R_0 = %lf\n",
                        __func__, solution->r_0);
#endif

        prec_t p = {
                .l = lower_bound_search(&eq, table, solution->r_0),
                .h = upper_bound_search(&eq, table, solution->r_0)
        };
#ifdef DEBUG
        printf("%s: solution: -R_0 = %lf\n", 
                        __func__, p.l);
        printf("%s: solution: +R_0 = %lf\n", 
                        __func__, p.h);
#endif

        get_errors(solution, table);
        iteration_storage_t *st = iteration_storage_create(table, solution);
        dump_all(solution, &p, st);

        // TODO: function's interface need improve
        if (cfg->mode == ITERATE_MODE)
                printf("_______%s_______\n", "First_stage_complete");
        cfg->filter = ERR_FILTER;
        cfg->h = get_limit_by_eps(table->size);
        cfg->l = sqrt(solution->sq / (solution->size + solution->s.size + 1));
#ifdef DEBUG
        printf("%s: kappa(%d) = %0.7lf\n", 
                        __func__,
                        table->size,
                        get_limit_by_eps(table->size));
#endif

        unsigned int old_size = table->size;
        while (cfg->mode == ITERATE_MODE) {
#ifdef DEBUG
                printf("%s: iter %d", __func__, j++);
#endif
                solution = __get_solution_iterate(table, &eq);
                cfg->h = get_limit_by_eps(table->size);
                cfg->l = sqrt(solution->sq / (solution->size + solution->s.size + 1));
                if (table->size == old_size)
                        break;
                old_size = table->size;
        }

        if (cfg->mode == ITERATE_MODE) {
                dump_objects_xyz(table, table->size);
                dump_table(table);
                dump_all(solution, &p, st);
        }
}

