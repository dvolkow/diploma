#ifndef CORE_H  
#define CORE_H  1

#define MAX_ORDER_SOLUTION      10

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

#endif // CORE_H
