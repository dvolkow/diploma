#ifndef IO_H
#define IO_H  1

#include "types.h"
#include "math.h"

#define INPUT_TABLE_FILE_NAME   \
        "apogee-rc.txt"

#define INPUT_FILE_PATH         \
        "./data/"

#define PRINT_IO_OPEN_ERROR(name)       \
        printf("%s: [io-error]: fail to open %s!\n",    \
                        __func__, name)


apogee_rc_table_t *read_table(const char *input_file_name);
void output_result(void);
unsigned int countlines(const char *filename);

opt_t *read_solution(const char *input_file_name);

#endif // IO_H

