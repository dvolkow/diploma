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


typedef enum {
        COUNTER,
        DISTANCE,
} averages_mode_t;

void dump_result(opt_t *, apogee_rc_table_t *, prec_t *);

#endif // GRAPH_H
