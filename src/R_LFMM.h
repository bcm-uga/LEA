/**
 * @file R_LFMM.h
 *
 * @brief C wrapper for LFMM command line program
 */

#ifndef R_LFMM_H
#define R_LFMM_H

void R_LFMM(char **R_input_file, char **R_output_file, char **R_cov_file,
            int *R_n, int *R_L, int *R_D, int *R_nd,
            int *R_K, int *R_Niter, int *R_burn,
            int *R_num_CPU, long long *R_seed, int *R_missing_data, int *R_all,
            double *R_dic, double *R_dev, double *noise_epsilon,
            double *b_epsilon, int *init);

#endif                          // R_LFMM_H
