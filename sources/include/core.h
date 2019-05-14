#ifndef CORE_H
#define CORE_H  1

#include <stdbool.h>
#include "filter.h"

#define LOWER_BOUND_R0          6.0
#define UPPER_BOUND_R0          12.0

#define W_SUN_START             7.0

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
        L_MODE,
        B_MODE,
} eq_mode_t;

typedef enum {
        SIMPLE_MODE
      , ITERATE_MODE
      , GENERATION_MODE
} g_mode_t;

typedef enum {
        VR_PART_MODE
      , B_PART_MODE
      , L_PART_MODE
      , UNI_MODE
      , UNINAT_MODE
      , INVALID_MODE
} solution_mode_t;


typedef struct {
#define DEFAULT_ORD             1
        unsigned int ord;
#define DEFAULT_INF_NAME        "apogee_rc.txt"
        char *input_file_name;
#define DEFAULT_DUMP_FILE_NAME  "dump_table.txt"
        char *dump_file_name;
#define DEFAULT_L               0
        double l;
#define DEFAULT_H               0
        double h;
        double n_err;
#define DEFAULT_MKSIZE          10
        unsigned int mksize;
#define DEFAULT_FILTER          BAD_FILTER
        filter_mode_t filter;
#define DEFAULT_MODE            SIMPLE_MODE
        g_mode_t mode;
#define DEFAULT_BOLTER          0
        int bolter;
#define DEFAULT_PROFILE         0
        int draw_profile;
	solution_mode_t solution_mode;
#define DEFAULT_SOL_MODE        INVALID_MODE
} parser_t;
#define GET_MODE(p_parser)		\
        ((p_parser)->mode)
#define GET_SOLUTION_MODE(p_parser)	\
        ((p_parser)->solution_mode)

void parse_args(int, char **);
bool parser_t_is_valid(const parser_t *cfg);
parser_t *get_parser(void);

int parser_init(void);
void parser_exit(void);


#define PREC_COMPARE_ITER                       1e-5
#define ITER_CONDITION(w_1, w_2, r_1, r_2)              \
        (fabs((w_1) - (w_2)) > PREC_COMPARE_ITER ||     \
         fabs((r_1) - (r_2)) > PREC_COMPARE_ITER)



#endif // CORE_H
