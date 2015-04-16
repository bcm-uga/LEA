/**
 * @file R_pca.h
 *
 * @brief C wrapper for pca command line program
 */

#ifndef R_pca_H
#define R_pca_H

void R_pca(char **R_input_file, char **R_output_eva_file,
           char **R_output_eve_file, char **R_output_sdev_file,
           char **R_output_x_file, int *n, int *L, int *K, int *c, int *s);

#endif                          // R_pca_H
