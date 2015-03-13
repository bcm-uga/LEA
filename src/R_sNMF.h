/**
 * @file R_sNMF.h
 *
 * @brief C wrapper for sNMF command line program
 */

#ifndef R_sNMF_H
#define R_sNMF_H

void R_sNMF(char **R_genotype_file, int *R_K, double *R_alpha, double *R_tol,
            double *R_percentage, int *R_iteration, long long *R_seed,
            int *R_ploidy, int *R_num_proc, char **R_input_file_Q,
            char **R_output_file_Q, char **R_output_file_G, int *I,
            double *all_ce, double *masked_ce, int *n, int *L);

#endif                          // R_sNMF_H
