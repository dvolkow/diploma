#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "core.h"
#include "mem.h"
#include "debug.h"
#include "filter.h"

#define NEXT_ARG(args, argv)    (argv++, args--)
#define CHECK_ARGS(args)        ((args) - 1 > 0)

#define __VERSION               "0.1"

static parser_t *__cfg;

parser_t *get_parser(void)
{
        return __cfg;
}

int parser_init(void)
{
        __cfg = dv_alloc(sizeof(parser_t));
        // TODO: default initializer?
        __cfg->ord              = DEFAULT_ORD;
        __cfg->input_file_name  = DEFAULT_INF_NAME;
        __cfg->dump_file_name  =  DEFAULT_DUMP_FILE_NAME;
        __cfg->l                = DEFAULT_L;
        __cfg->h                = DEFAULT_H;
        __cfg->filter           = DEFAULT_FILTER;
        __cfg->mode             = DEFAULT_MODE;
        __cfg->mksize           = DEFAULT_MKSIZE;
        __cfg->bolter           = DEFAULT_BOLTER;

        return 0;
}


static void show_usage()
{
        printf("\e[1mSYNOPSIS\e[0m\n");
        printf("                \e[1mdiploma\e[0m -f IN_FILE --ord ORD [\e[4m--filter\e[0m] \e[4mFILTER L H\e[0m \n");
        printf("For more information and samples type 'diploma -h'.\n\n");
        printf("\e[1mAUTHOR\e[0m\n");
        printf("               Written by Daniel Wolkow (2018).\n");
        printf("\e[1mVERSION\e[0m\n");
        printf("               Your copy is %s, ", __VERSION);
        printf("newest version by address \n");
        printf("               \e[4mhttps://github.com/dvolkow/diploma\e[0m\n");
}

bool parser_t_is_valid(const parser_t *cfg)
{
        return cfg->input_file_name != NULL &&
               cfg->ord > 0;
}

static inline int matches(const char *arg,
                          const char *argv)
{
        return !strcmp(argv, arg);
}


static void filter_assigner(const char *argv)
{
        if (matches("L", argv)) {
                __cfg->filter = L_FILTER;
        } else if (matches("B", argv)) {
                __cfg->filter = B_FILTER;
        } else if (matches("ERR", argv)) {
                __cfg->filter = ERR_FILTER;
        } else {
                __cfg->filter = BAD_FILTER;
        }
}

static void mode_assigner(const char *argv)
{
        if (matches("I", argv)) 
                __cfg->mode = ITERATE_MODE;
        if (matches("G", argv))
                __cfg->mode = GENERATION_MODE;
}

void parse_args(int argc, 
                char *argv[])
{
        parser_t *res = get_parser();
        double readbuff;

        while (CHECK_ARGS(argc)) {
                NEXT_ARG(argc, argv);
                if (matches("-f", *argv) || matches("--file", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                res->input_file_name = *argv;
                        } else {
                                goto usage_ret;
                        }
                } else if (matches("--ord", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                res->ord = (unsigned int)atoi(*argv);
                        } else {
                                goto usage_ret;
                        }
                } else if (matches("--mksize", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                res->mksize = (unsigned int)atoi(*argv);
                        } else {
                                goto usage_ret;
                        }
                } else if (matches("-d", *argv) || matches("--dump", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                res->dump_file_name = *argv;
                        } else {
                                goto usage_ret;
                        }
                } else if (matches("-e", *argv) || matches("--err", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                sscanf(*argv, "%lf", &readbuff);
                                res->n_err = readbuff;
                        } else {
                                goto usage_ret;
                        }
                } else if (matches("-m", *argv) || matches("--mode", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                mode_assigner(*argv);
                        } else {
                                goto usage_ret;
                        }
                } else if (matches("-b", *argv) || matches("--bolter", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                res->bolter = atoi(*argv);
                        } else {
                                goto usage_ret;
                        }
                } else if (matches("--filter", *argv)) {
                        if (CHECK_ARGS(argc)) {
                                NEXT_ARG(argc, argv);
                                filter_assigner(*argv);
                                if (CHECK_ARGS(argc)) {
                                        NEXT_ARG(argc, argv);
                                        sscanf(*argv, "%lf", &readbuff);
                                        res->l = readbuff;
                                }
                                if (CHECK_ARGS(argc)) {
                                        NEXT_ARG(argc, argv);
                                        sscanf(*argv, "%lf", &readbuff);
                                        res->h = readbuff;
                                } 
                        } else {
                                goto usage_ret;
                        }
                } 
        }

        return;
usage_ret:
        show_usage();
}






void parser_exit() {}
