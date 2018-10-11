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

void dump_result(opt_t *, apogee_rc_table_t *, prec_t *);
void dump_averages(iteration_storage_t *, opt_t *, averages_mode_t);

#define DEFAULT_BACKGROUND_COUNT        42

void dump_background(const iteration_storage_t *, 
                     const opt_t *, 
                     const int);
#endif // GRAPH_H
