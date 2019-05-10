#include <assert.h>
#include "db.h"
#include "types.h"
#include "generators.h"
#include "mem.h"


static apogee_rc_table_t **db;
static unsigned int size_db;


int db_init(void)
{
        db = dv_alloc(sizeof(apogee_rc_table_t *) * DB_SIZE);
        size_db = 0;

        return 0;
}


void db_add(apogee_rc_table_t *table)
{
        if (size_db < DB_SIZE)
                db[size_db++] = table;
        else
                printf("%s: warning %d\n", __func__, __LINE__);

        table->size = 0;
}

unsigned int db_size(void)
{
        return size_db;
}

void db_add_to_table(apogee_rc_t *line, unsigned int type)
{
        assert(type < db_size());
        unsigned int cur_top = db[type]->size;
        assert(cur_top < PREALLOC_TABLE_SIZE);
        db[type]->data[cur_top] = *line;
        db[type]->size++;
}


apogee_rc_table_t *db_get(unsigned int type)
{
        assert(type < db_size());
        return db[type];
}

void db_exit(void)
{
}
