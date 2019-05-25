#ifndef IO_H
#define IO_H  1

#include "types.h"
#include "math.h"

#define INPUT_TABLE_FILE_NAME   \
        "apogee-rc.txt"

#define INPUT_FILE_PATH         \
        "./data/"

#define PRINT_IO_OPEN_ERROR(name)       \
        printf("%s: [IOERROR]: fail to open %s!\n",    \
                        __func__, name)

#define CHECK_FILE_AND_RET(fdesc, name)                 \
        do {                                            \
                if ((fdesc) == NULL) {                  \
                        PRINT_IO_OPEN_ERROR(name);      \
                        return;                         \
                }                                       \
        } while (0)

#define K_A     1.044 // Calibration coeff.


apogee_rc_table_t *read_table(const char *input_file_name);
void dump_table(const apogee_rc_table_t *);
unsigned int countlines(const char *filename);

opt_t *read_solution(const char *input_file_name);

#endif // IO_H

