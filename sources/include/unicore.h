#ifndef UNICORE_H
#define UNICORE_H  1

#include "types.h"
#include "opt.h"


#define VR_PART         0
#define L_PART          1
#define B_PART          2
#define TOTAL_QTY       3

double core_vr_get_beta_n(const apogee_rc_t *, beta_ord_t);
double core_vr_get_alpha_n(const apogee_rc_t *,
                           const unsigned int ,
                           const double );
void uni_g_sd_init(const double *);
void vr_b_iterations(apogee_rc_table_t *);
opt_t *united_entry(apogee_rc_table_t *);
opt_t *united_with_nature_errs_entry(apogee_rc_table_t *);

void precalc_errors_uni(apogee_rc_table_t *table,
                        const double limit);
void filter_get_and_apply(apogee_rc_table_t *table);
void get_solution(apogee_rc_table_t *table);
void get_partial_vr_solution(apogee_rc_table_t *table);

#endif // UNICORE_H  1
