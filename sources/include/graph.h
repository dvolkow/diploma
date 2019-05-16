#ifndef GRAPH_H
#define GRAPH_H         1

#include "types.h"
#include "io.h"
#include "opt.h"

#define EDGE_SPARCE_DIVISOR     10
#define OMEGA_SUN               30.50
#define ROTC_STEP_R             0.1
#define ROTC_LOWER_BOUND        2
#define ROTC_UPPER_BOUND        15

#define MAX_PRINTF_COLS		12

#define RC_OUT_FILE_NAME        \
        "rotc.txt"

#define DF_OUT_FILE_NAME        \
        "objs.txt"

#define SUN_POINT_FILE_NAME     \
        "sun.txt"

#define AVERAGE_R_FILE_NAME     \
        "averages.txt"

#define R0THETA0_FILE_NAME      \
        "R0Theta0.txt"

#define R0THETA0_MAIN_FILE_NAME \
        "R0Theta0_main.txt"

#define AVERAGE_COUNT_BASE              500
#define AVERAGE_COUNT_EDGE              (AVERAGE_COUNT_BASE / EDGE_SPARCE_DIVISOR)
#define AVERAGE_COUNT_EDGE_L_BOUND      5
#define AVERAGE_COUNT_EDGE_R_BOUND      12

#define IS_LAST_STEP(c, s)                              \
        (((c) - (s)) <= 0)

#define R_INTO_MIDDLE(v)                                \
        ((v < AVERAGE_COUNT_EDGE_R_BOUND) &&            \
                        (v > AVERAGE_COUNT_EDGE_L_BOUND))



typedef enum {
        COUNTER,
        DISTANCE,
} averages_mode_t;


typedef struct {
        double r;
        double theta;
        // bounds from Monte-Karlo:
        double theta_min;
        double theta_max;
} rot_curve_t;

#define DEFAULT_BACKGROUND_COUNT        127
void dump_rand_test(const double *,
                    const dsize_t);

void dump_objects_xyz(const apogee_rc_table_t *,
                      const dsize_t,
                      const char *);
void dump_result(const opt_t *opt);

void dump_core_l_solution(const opt_t *);
void dump_core_b_solution(const opt_t *);
void dump_core_vr_solution(const opt_t *);
void dump_united_solution(const opt_t *);
void dump_united_solution_points(const opt_t *);

void dump_uni_rotation_curve(const rot_curve_t *, const unsigned int);
void dump_uni_rotation_objs(const apogee_rc_table_t *,
                            const opt_t *);

void dump_rotation_curve_vr(const opt_t *solution);
void dump_rotation_curve_b(const opt_t *solution);
void dump_rotation_curve_l(const opt_t *solution);

double get_point_by_uni_solution(const opt_t *,
                                 const double);
double get_point_by_vr_solution(const opt_t *,
                                const double);
double get_point_by_b_solution(const opt_t *,
                                const double);
double get_point_by_l_solution(const opt_t *,
                                const double);

void dump_table_parameters(const apogee_rc_table_t *,
                           const opt_t *);

void dump_objects_theta_R(const apogee_rc_table_t *,
                          opt_t *,
                          unsigned int,
                          const char *);

double get_c_point_by_part_solution(const opt_t *,
                                    const double,
                                    unsigned int,
                                    const double,
                                    unsigned int);

void dump_part_rotation_curve(const opt_t *,
                              unsigned int,
                              const char *,
                              const double);

void dump_vr_solution(const opt_t *);
void dump_R0_theta_ellips(const opt_t **,
                          unsigned int,
                          const opt_t *);

char *name_for_obj(const unsigned int i,
                   const unsigned int n,
                   const char *prefix);

void dump_profile(linear_equation_t *,
                  apogee_rc_table_t *,
                  opt_params_t *,
                  const char *);

void multiply_dump_unfriendly_result(const opt_t **,
				     const unsigned int,
				     const char *);
void partial_dump_unfriendly_result(const opt_t *,
				    unsigned int,
				    const char *);
#endif // GRAPH_H
