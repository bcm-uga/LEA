/**
 * @addtogroup register_pca
 * @ingroup pca
 * @{
 * @file register_pca.h
 *
 * @brief functions to analyze the command-line
 */

#ifndef REGISTER_PCA_H
#define REGISTER_PCA_H

/**
 * analyse command line set of parameters and set the parameters
 * 
 * @param argc  	the number of arguments
 * @param argv  	the set of arguments
 * @param input		the input file
 * @param output_eva	eigenvalues file	
 * @param output_eve	eigenvectors file	
 * @param output_sdev	standard deviation file	
 * @param output_x	projection file	
 * @param K     	the number of clusters
 * @param c     	boolean: data centered
 * @param s     	boolean: data scaled
 */
void analyse_param_pca( int argc, char *argv[], char *input, char* output_eva,
                        char* output_eve, char *output_sdev, char *output_x, 
			int *K, int *c, int *s);

#endif // REGISTER_PCA_H

/** @} */
