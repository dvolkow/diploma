#ifndef IO_H  
#define IO_H  1

#include "types.h"

#define INPUT_TABLE_FILE_NAME   \
        "apogee-rc.txt"

#define INPUT_FILE_PATH         \
        "./data/"

apogee_rc_table_t *read_table(const char *input_file_name);
void output_result(void);
unsigned int countlines(char *filename);


#endif // IO_H

