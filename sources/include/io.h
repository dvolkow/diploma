#ifndef IO_H  
#define IO_H  1

#include "core.h"

#define INPUT_TABLE_FILE_NAME   \
        "apogee-rc.txt"

#define INPUT_FILE_PATH         \
        "./data/"

/**
 * TODO: description for fields
 *
 */
typedef struct {
        double l;
        double b;
        double v_helio;
        double dist;
        double pm_ra;
        double pm_dec;
} apogee_rc_t;


typedef struct {
        double err[MAX_ORDER_SOLUTION]; 
        double r; // R 
        double theta;
} iteration_storage_t;

apogee_rc_t *read_table(void);
void output_result(void);
unsigned int countlines(char *filename);


#endif // IO_H

