#ifndef DB_H
#define DB_H    1
#include "types.h"

#define DB_SIZE         16

enum {
        ERROR_LIMITED = 0,
        MATCH_LIMITED,
};

int db_init(void);
void db_exit(void);
void db_add(apogee_rc_table_t *);
unsigned int db_size(void);
apogee_rc_table_t *db_get(unsigned int);

#endif // DB_H
