#ifndef FILTER_H
#define FILTER_H         1

typedef enum {
        L_FILTER
      , B_FILTER
      , ERR_FILTER
      , UNITED_FILTER
      , MATCH_FILTER

      // must be last filter:
      , BAD_FILTER 
} filter_mode_t;

/**
 * @l: lower bound for parameter
 * @h: upper bound for parameter
 * @f: oracle for line of table that 
 *     const void* parameter
 */
typedef struct {
        double l;
        double h;
        bool (*f)(const void *, const double, const double);
} filter_t;


#endif // FILTER_H

