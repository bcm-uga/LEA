/**
 * @file R_crossEntropy.h
 *
 * @brief C wrapper for crossEntropy command line program
 */

#ifndef R_crossEntropy_H
#define R_crossEntropy_H

void R_crossEntropy(char **R_genotype_file, char **R_missing_data_file,
                    char **R_Q_file, char **R_F_file, int *R_K, int *R_ploidy,
                    double *all_ce, double *masked_ce);

#endif                          // R_crossEntropy_H
