#ifndef CORE_H  
#define CORE_H  1

#include <stdbool.h>
#include "filter.h"

#define LOWER_BOUND_R0          3.0
#define UPPER_BOUND_R0          17.0

#define SEARCH_PRECISION        0.001
#define STEP_DIVISOR            4

typedef enum {
        FIRST = 0,
        SECOND,
        THIRD,

        //  MUST BE LAST:
        BETA_QTY,
} beta_ord_t;


typedef enum {
        VR_MODE,
        PM_RA_MODE,
        PM_DEC_MODE,
} eq_mode_t;

typedef enum {
        SIMPLE_MODE
      , ITERATE_MODE
      , GENERATION_MODE
} g_mode_t;

typedef struct {
#define DEFAULT_ORD             1
        int ord;
#define DEFAULT_INF_NAME        "apogee_rc.txt"
        char *input_file_name;
#define DEFAULT_DUMP_FILE_NAME  "dump_table.txt"
        char *dump_file_name;
#define DEFAULT_L               0
        double l;
#define DEFAULT_H               0
        double h;
#define DEFAULT_FILTER          BAD_FILTER
        filter_mode_t filter;
#define DEFAULT_MODE            SIMPLE_MODE
        g_mode_t mode;
} parser_t;
#define GET_MODE(p_parser)      \
        ((p_parser)->mode)

void parse_args(int, const char **);
bool parser_t_is_valid(const parser_t *cfg);
parser_t *get_parser(void);

void get_solution(void);
int parser_init(void);
void parser_exit(void);
#endif // CORE_H
