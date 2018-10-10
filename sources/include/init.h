#ifndef INIT_H
#define INIT_H          1

typedef struct {
        char *name;
        int (*init)(void);
        void (*exit)(void);
} init_t;

#endif // INIT_H
