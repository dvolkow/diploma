#ifndef CORE_H  
#define CORE_H  1

#include <stdbool.h>

#define MAX_ORDER_SOLUTION      10

#define LOWER_BOUND_R0          3.0
#define UPPER_BOUND_R0          17.0

#define SEARCH_PRECISION        0.001
#define STEP_DIVISOR            2

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
        NAME_PROGRAM_ARG = 0,
        ORD_ARG,
        INPUT_FILE_ARG,

        // MUST BE LAST:
        CMD_ARGS_QTY
} cmd_line_arg_t;


typedef struct {
#define DEFAULT_ORD             1
        int ord;
#define DEFAULT_INF_NAME        "apogee_rc.txt"
        char *input_file_name;
#define DEFAULT_L               0
        double l;
#define DEFAULT_H               0
        double h;
#define DEFAULT_FILTER          NULL
        char *filter;
} parser_t;

parser_t *parse_args(int, const char **);
bool parser_t_is_valid(const parser_t *cfg);

#endif // CORE_H
