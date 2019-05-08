#ifndef GRAPH_H
#define GRAPH_H         1

#include "types.h"
#include "io.h"

#define EDGE_SPARCE_DIVISOR     10
#define OMEGA_SUN               30.50
#define ROTC_STEP_R             0.1
#define ROTC_LOWER_BOUND        2
#define ROTC_UPPER_BOUND        15

#define RC_OUT_FILE_NAME        \
        "rotc.txt"

#define DF_OUT_FILE_NAME        \
        "objs.txt"

#define SUN_POINT_FILE_NAME     \
        "sun.txt"

#define AVERAGE_R_FILE_NAME     \
        "averages.txt"

#define AVERAGE_COUNT_BASE              500
#define AVERAGE_COUNT_EDGE              (AVERAGE_COUNT_BASE / EDGE_SPARCE_DIVISOR)
#define AVERAGE_COUNT_EDGE_L_BOUND      5
#define AVERAGE_COUNT_EDGE_R_BOUND      12

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
void dump_rand_test(const double *array, 
                    const dsize_t size);

void dump_table(const apogee_rc_table_t *table);
void dump_objects_xyz(const apogee_rc_table_t *table, const dsize_t size);
void dump_all(opt_t *solution, prec_t *p, iteration_storage_t *st);

void dump_core_l_solution(const opt_t *);
void dump_core_b_solution(const opt_t *);
void dump_core_vr_solution(const opt_t *);
void dump_united_solution(const opt_t *);

void dump_uni_rotation_curve(const rot_curve_t *, const unsigned int);
void dump_uni_rotation_objs(const apogee_rc_table_t *,
                            const opt_t *);

double get_point_by_uni_solution(const opt_t *, const double);
void dump_table_parameters(const apogee_rc_table_t *,
                           const opt_t *);
#endif // GRAPH_H
