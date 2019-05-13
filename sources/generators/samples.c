#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <assert.h>

#include "math.h"
#include "mem.h"
#include "core.h"
#include "opt.h"
#include "generators.h"
#include "graph.h"
#include "unicore.h"
#include "trigonometry.h"

static gsl_rng *g_gauss_rg_p;


double *gen_vector_by_mean_and_sd(const gsl_rng *r, 
                                  const double mean, 
                                  const double sigma,
                                  const unsigned int size)
{
        unsigned int i;
        double *res = dv_alloc(size);
        for (i = 0; i < size; ++i) {
                res[i] = mean + gsl_ran_gaussian(r, sigma);
        }

        return res;
}

static inline double __gen_gauss_d(const gsl_rng *r, 
                                   const double mean, 
                                   const double sigma)
{
        return mean + gsl_ran_gaussian(r, sigma);
}

static inline double __gen_flat_d(const gsl_rng *r,
                                 const double l,
                                 const double h)
{
        return gsl_ran_flat(r, l, h);
}

static double __gen_dist_d(const gsl_rng *r,
                           const double mean, 
                           const double sd)
{
        double res;
        do {
                res = fabs(__gen_gauss_d(r, mean, sd));
        } 
        while (res > 12);
        return res;
}

static double __gen_b_d(const gsl_rng *r,
                           const double mean,
                           const double sd)
{
#define CUT_BY_B        80
        double res;
        do {
                res = __gen_gauss_d(r, mean, sd);
        }
        while (fabs(res) >= CUT_BY_B);

        return res;
}

void fill_table_by_vr_solution(const opt_t *solution,
                               const apogee_rc_table_t *src_table,
                               apogee_rc_table_t *dst_table)
{
        gsl_rng *rng = dv_rand_acquire(MT);
        assert(rng != NULL);

        dst_table->r_0 = solution->r_0;
        dst_table->size = src_table->size;

        const double sd = sqrt(src_table->sigma[VR_PART]);
        unsigned int i;
        for (i = 0; i < src_table->size; ++i) {
                dst_table->data[i].l = src_table->data[i].l;
                dst_table->data[i].b = src_table->data[i].b;
                dst_table->data[i].cos_l = src_table->data[i].cos_l;
                dst_table->data[i].cos_b = src_table->data[i].cos_b;
                dst_table->data[i].sin_b = src_table->data[i].sin_b;
                dst_table->data[i].sin_l = src_table->data[i].sin_l;
                dst_table->data[i].dist = src_table->data[i].dist;

                dst_table->data[i].v_helio = __gen_gauss_d(rng,
                                                get_mod_vr(solution, &src_table->data[i]),
                                                  sd);
        }
        dv_rand_release(rng);
}

void fill_table_by_b_solution(const opt_t *solution,
                               const apogee_rc_table_t *src_table,
                               apogee_rc_table_t *dst_table)
{
        gsl_rng *rng = dv_rand_acquire(MT);
        assert(rng != NULL);

        dst_table->r_0 = solution->r_0;
        dst_table->size = src_table->size;

        const double sd = sqrt(src_table->sigma[B_PART]);
        unsigned int i;
        for (i = 0; i < src_table->size; ++i) {
                dst_table->data[i].l = src_table->data[i].l;
                dst_table->data[i].b = src_table->data[i].b;
                dst_table->data[i].cos_l = src_table->data[i].cos_l;
                dst_table->data[i].cos_b = src_table->data[i].cos_b;
                dst_table->data[i].sin_b = src_table->data[i].sin_b;
                dst_table->data[i].sin_l = src_table->data[i].sin_l;
                dst_table->data[i].dist = src_table->data[i].dist;

                dst_table->data[i].pm_b = __gen_gauss_d(rng,
                                                get_mod_b(solution, &src_table->data[i]),
                                                  sd);
        }
        dv_rand_release(rng);
}

void fill_table_by_l_solution(const opt_t *solution,
                               const apogee_rc_table_t *src_table,
                               apogee_rc_table_t *dst_table)
{
        gsl_rng *rng = dv_rand_acquire(MT);
        assert(rng != NULL);

        dst_table->r_0 = solution->r_0;
        dst_table->size = src_table->size;

        const double sd = sqrt(src_table->sigma[L_PART]);
        unsigned int i;
        for (i = 0; i < src_table->size; ++i) {
                dst_table->data[i].l = src_table->data[i].l;
                dst_table->data[i].b = src_table->data[i].b;
                dst_table->data[i].cos_l = src_table->data[i].cos_l;
                dst_table->data[i].cos_b = src_table->data[i].cos_b;
                dst_table->data[i].sin_b = src_table->data[i].sin_b;
                dst_table->data[i].sin_l = src_table->data[i].sin_l;
                dst_table->data[i].dist = src_table->data[i].dist;

                dst_table->data[i].pm_l = __gen_gauss_d(rng,
                                                get_mod_l(solution, &src_table->data[i]),
                                                  sd);
        }
        dv_rand_release(rng);
}



void fill_table_by_uni_solution(const opt_t *solution,
                                const apogee_rc_table_t *src_table,
                                apogee_rc_table_t *dst_table)
{
        gsl_rng *rng = dv_rand_acquire(MT);
        assert(rng != NULL);

        const unsigned int ssize = solution->size;
        dst_table->r_0 = solution->r_0;
        dst_table->size = ssize;

        unsigned int i;
        double sd[TOTAL_QTY];
        for (i = 0; i < TOTAL_QTY; ++i) {
                sd[i] = sqrt(src_table->sigma[i]);
#ifdef DEBUG_GEN
                printf("%s: %u = %lf\n", __func__, i, sd[i]);
#endif
        }

        for (i = 0; i < ssize; ++i) {
                //dst_table->data[i].id = i;
                dst_table->data[i].l = src_table->data[i].l;
                dst_table->data[i].b = src_table->data[i].b;
                dst_table->data[i].cos_l = src_table->data[i].cos_l;
                dst_table->data[i].cos_b = src_table->data[i].cos_b;
                dst_table->data[i].sin_b = src_table->data[i].sin_b;
                dst_table->data[i].sin_l = src_table->data[i].sin_l;
                dst_table->data[i].dist = src_table->data[i].dist;
                dst_table->data[i].pm_match = src_table->data[i].pm_match;

//              TODO: Variable for solutions:
//              dst_table->data[i].pm_b_err = src_table->data[i].pm_b_err;
//              dst_table->data[i].pm_l_err = src_table->data[i].pm_l_err;

                dst_table->data[i].v_helio = __gen_gauss_d(rng,
                                                       src_table->data[i].v_helio,
                                                       sd[VR_PART]);
#ifdef DEBUG_GEN
                if (i < 10)
                printf("%s[VR]: origin %lf, new %lf\n", __func__, 
                                src_table->data[i].v_helio, 
                                dst_table->data[i].v_helio);
#endif
                dst_table->data[i].pm_l = __gen_gauss_d(rng,
                                                       src_table->data[i].pm_l,
                                                       sd[L_PART]);
#ifdef DEBUG_GEN
                if (i < 10)
                printf("%s[PM_L]: origin %lf, new %lf\n", __func__, 
                                src_table->data[i].pm_l, 
                                dst_table->data[i].pm_l);
#endif
                dst_table->data[i].pm_b = __gen_gauss_d(rng,
                                                       src_table->data[i].pm_b,
                                                       sd[B_PART]);
#ifdef DEBUG_GEN
                if (i < 10)
                printf("%s[PM_B]: origin %lf, new %lf\n", __func__,
                                src_table->data[i].pm_b,
                                dst_table->data[i].pm_b);
#endif
        }

        dv_rand_release(rng);
        vr_b_iterations(dst_table);
}



opt_t *monte_carlo_entry(const opt_t *solution,
                         const apogee_rc_table_t *data,
                         const mk_params_t *params)
{
        const unsigned int count = params->count;
        const opt_t **results = dv_alloc(sizeof(opt_t *) * count);
        apogee_rc_table_t *tmp_table = dv_alloc(sizeof(apogee_rc_table_t));
        const unsigned int ssize = data->size;
        const unsigned int n = solution->s.size;
        tmp_table->data = dv_alloc(sizeof(apogee_rc_t) * ssize);
        tmp_table->r_0 = solution->r_0;

        opt_t main_res = {
                .s = { 0 },
                .sq = solution->sq,
                .size = data->size
        };
        main_res.s.data = dv_alloc(sizeof(double) * n);
        main_res.s.size = n;
        main_res.bounds = dv_alloc(sizeof(prec_t) * n);

        const unsigned int fits_count = (unsigned int)((ROTC_UPPER_BOUND - ROTC_LOWER_BOUND) / ROTC_STEP_R);
        rot_curve_t *curve = dv_alloc(sizeof(rot_curve_t) * fits_count);
        double **__curve = dv_alloc(sizeof(double *) * count);
        double r = ROTC_LOWER_BOUND;
        unsigned int i, j;

#ifdef DEBUG_GEN
        double *test_vr = dv_alloc(sizeof(double) * count);
        double *test_l = dv_alloc(sizeof(double) * count);
        double *test_b = dv_alloc(sizeof(double) * count);
#endif

        i = 0;
        while (r < ROTC_UPPER_BOUND) {
                curve[i].theta = params->f_point_by_solution(solution, r);
                curve[i].theta_max = curve[i].theta;
                curve[i].theta_min = curve[i].theta;
                curve[i].r = r;
                ++i;
                r += ROTC_STEP_R;
        }

        for (i = 0; i < count; ++i) {
                params->f_table_by_solution(solution, data, tmp_table);

#ifdef DEBUG_GEN
                test_vr[i] = tmp_table->data[15].v_helio;
                test_l[i] = tmp_table->data[15].pm_l;
                test_b[i] = tmp_table->data[15].pm_b;
#endif

                results[i] = params->f_entry(tmp_table);
                printf("%s: [%d/%d] completed\n", __func__, i + 1, count);
                //dump_united_solution_points(results[i]);
                __curve[i] = (double *)dv_alloc(sizeof(double) * fits_count);
                for (j = 0; j < fits_count; ++j) {
                        double theta = params->f_point_by_solution(results[i], curve[j].r);
                        __curve[i][j] = theta;
                }
        }

#ifdef DEBUG_GEN
        printf("%s: test_vr %lf pm %lf, test_l %lf pm %lf, test_b %lf pm %lf\n", __func__,
                        get_mean(test_vr, count), get_sd(test_vr, count),
                        get_mean(test_l, count), get_sd(test_l, count),
                        get_mean(test_b, count), get_sd(test_b, count));
#endif

        //dump_uni_rotation_objs(data, solution);

        double *tmp_line = dv_alloc(sizeof(double) * count);
        for (j = 0; j < n; ++j) {
                for (i = 0; i < count; ++i) {
                        tmp_line[i] = results[i]->s.data[j];
                }
                main_res.bounds[j].l = get_sd(tmp_line, count);
                main_res.s.data[j] = get_mean(tmp_line, count);
        }

        for (i = 0; i < fits_count; ++i) {
                for (j = 0; j < count; ++j)
                        tmp_line[j] = __curve[j][i];
                double sd_curve = get_sd(tmp_line, count);
                curve[i].theta_max = curve[i].theta + sd_curve;
                curve[i].theta_min = curve[i].theta - sd_curve;
        }

        dump_uni_rotation_curve(curve, fits_count);

        for (j = 0; j < count; ++j) {
                tmp_line[j] = results[j]->r_0;
        }

        main_res.r_0 = get_mean(tmp_line, count);
        main_res.dr_0 = get_sd(tmp_line, count);
        for (j = 0; j < count; ++j) {
                tmp_line[j] = results[j]->sq;
        }

        dump_R0_theta_ellips(results, count, solution);

        main_res.sq = get_mean(tmp_line, count);
        opt_t *ret = dv_alloc(sizeof(opt_t));
        *ret = main_res;
        /**
         * Test:
        double *tmp_nums = dv_alloc(sizeof(double) * count);
        gsl_rng *rng = dv_rand_acquire(MT);
        for (i = 0; i < count; i++) {
                tmp_nums[i] = __gen_gauss_d(rng, 10, 32);
        }
        printf("%s: mean %lf, sd %lf\n", __func__, get_mean(tmp_nums, count),
                                                   get_sd(tmp_nums, count));
        dv_rand_release(rng);
         */
        return ret;
}


apogee_rc_table_t *gen_table_by_solution(const opt_t *solution)
{
        apogee_rc_table_t *table = dv_alloc(sizeof(apogee_rc_table_t));
        gsl_rng *rng = dv_rand_acquire(MT);
        assert(rng != NULL);


        const unsigned int ssize = solution->size;
        table->r_0 = solution->r_0;
        table->size = ssize;
        table->data = dv_alloc(sizeof(apogee_rc_t) * ssize);

        dsize_t i;
#define LLOW    0
#define LHIGH   360

#define BLOW    -90
#define BHIGH   90

#define RLOW    0.39
#define RHIGH   12.83    

        for (i = 0; i < ssize; ++i) {
                //table->data[i].id = i;
                table->data[i].l = deg_to_rad(__gen_flat_d(rng, LLOW, LHIGH));
                table->data[i].b = deg_to_rad(__gen_b_d(rng, 0, 15)); //__gen_flat_d(rng, BLOW, BHIGH);
                table->data[i].dist = __gen_dist_d(rng, 3, 1.5); //__gen_flat_d(rng, RLOW, RHIGH);
                table->data[i].v_helio = __gen_gauss_d(rng, get_mod_vr(solution, 
                                                                &table->data[i]), solution->sq);
        }

        dv_rand_release(rng);
        return table;
}

double *gen_vector_by_bounds_uni(const gsl_rng *r,
                                 const double l,
                                 const double h,
                                 const unsigned int size)
{
        unsigned int i;
        double *res = dv_alloc(size);

        for (i = 0; i < size; ++i) {
                res[i] = __gen_flat_d(r, l, h);
        }

        return res;
}

static gsl_rng *__dv_rand_acquire(const unsigned long int seed,
                                  gen_mode_t mode)
{
        gsl_rng *r;
        switch (mode) {
                case MT:
                default:
                        r = gsl_rng_alloc(gsl_rng_mt19937);
                        break;
        }
        gsl_rng_set(r, seed);

        return r;
}

gsl_rng *dv_rand_acquire(gen_mode_t mode)
{
#define get_random_double_by_memory()  \
        (*(double *)dv_mm_get_current_top())
#ifdef DEBUG
        printf("%s: get_random_double_by_memory = %lf\n",
                        __func__, get_random_double_by_memory());
#endif
        return __dv_rand_acquire((unsigned long int)(gsl_rng_uniform(g_gauss_rg_p) 
                                        * 1000), mode);
}

void dv_rand_release(gsl_rng *r)
{
        gsl_rng_free(r);
}


apogee_rc_table_t *create_apogee_rc_table_by_size(unsigned int size)
{
        apogee_rc_t *apogee_rc = dv_alloc(sizeof(apogee_rc_t) * size);
        apogee_rc_table_t *table = dv_alloc(sizeof(apogee_rc_table_t));
        table->data = apogee_rc;
        table->size = size;
        return table;
}

apogee_rc_table_t *generic_table(void)
{
        return create_apogee_rc_table_by_size(PREALLOC_TABLE_SIZE);
}


void generate(void)
{
        parser_t *cfg = get_parser();

        opt_t *solution = read_solution(cfg->input_file_name);
        if (solution == NULL) {
                printf("%s: failure when read solution\n", __func__);
                return;
        }

        apogee_rc_table_t *table = gen_table_by_solution(solution);
        dump_table(table);
        dump_objects_xyz(table, table->size, name_for_obj(__LINE__, 0, "xyz"));
}

int random_seed_init()
{
        g_gauss_rg_p = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(g_gauss_rg_p, dv_random_seed());
        return 0;
}

void random_seed_exit()
{
        dv_rand_release(g_gauss_rg_p);
}
