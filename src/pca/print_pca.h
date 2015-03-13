/**
 * @addtogroup print_pca
 * @ingroup pca
 * @{
 * @file print_pca.h
 *
 * @brief functions to print help and summary of the parameters
 */


#ifndef PRINT_PCA_H
#define PRINT_PCA_H

/**
 * print help
 */
void print_help_pca();


/**
 * print summary of the parameters
 *
 * @param N     the number of individuals
 * @param M     the number of loci
 * @param K     the number of latent factors
 * @param c	data centered (boolean)
 * @param s	data scaled (boolean)
 * @param input         genotype file
 * @param output_values eigenvalues file
 * @param output_vectors eigenvectors file
 * @param sdev_values 	standard deviation file
 * @param x_vectors	projection file
 */
void print_summary_pca( int N, int M, int K, int c, int s,
                        char *input, char *output_values,
                        char *output_vectors, char *sdev_value,
			char *x_vectors);    

#endif // PRINT_PCA_H

/** @} */
