#ifndef UNICORE_H
#define UNICORE_H  1

#include "types.h"
#include "opt.h"


#define VR_PART         0
#define L_PART          1
#define B_PART          2
#define TOTAL_QTY       3

double core_vr_get_beta_n(const apogee_rc_t *,
                          beta_ord_t);
double core_vr_get_alpha_n(const apogee_rc_t *,
                           const unsigned int ,
                           const double );
void uni_g_sd_init(const double *);
void vr_b_iterations(apogee_rc_table_t *);
opt_t *united_entry(apogee_rc_table_t *);
opt_t *united_with_nature_errs_entry(apogee_rc_table_t *);

void precalc_errors_uni(apogee_rc_table_t *,
                        const double);
void precalc_errors_uni_sigma_0(apogee_rc_table_t *,
                        const double);
void filter_get_and_apply(apogee_rc_table_t *);

void get_united_solution(apogee_rc_table_t *);
void get_united_sigma_0_solution(apogee_rc_table_t *);
void find_united_sigma_0_solution(apogee_rc_table_t *);

void get_partial_vr_solution(apogee_rc_table_t *);

void precalc_vsd_to_dump(apogee_rc_table_t *);


double get_v_generic_from_uni(const linear_eq_solve_t *,
                              const apogee_rc_t *,
                              const double,
                              const unsigned int);
void dump_united_solution_profile(apogee_rc_table_t *,
                                  unsigned int);
void dump_united_sigma0_solution_profile(apogee_rc_table_t *,
                                         unsigned int);
double sigma_for_k(const unsigned int,
                   const apogee_rc_t *,
                   const double);

void dump_united_solution_r0_bounds(apogee_rc_table_t *,
				    opt_t *);
#endif // UNICORE_H  1
