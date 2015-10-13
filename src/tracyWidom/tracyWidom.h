/**
 * @addtogroup tracyWidom
 * @ingroup tracyWidom
 * @{
 * @file tracyWidom.h
 *
 * @brief functions to perform tracy-widom test on a set of eigenvalues
 * 	  such as described by Patterson et al (2006) "Eigenanalysis".
 *	  This program is very similar with the one in Eigensoft.
 */

#ifndef TRACYWIDOM_H
#define TRACYWIDOM_H

/**
 * run tracy widom tests of eigenvalues
 *
 * @param input_file	input file name of eigenvalues
 * @param output_file   the output file with p-values
 */
void tracyWidom(char *input_file, char *output_file);

/**
 * write data into output file
 *
 * @param output_file	output file name
 * @param M		number of non-zero eigenvalues
 * @param values	eigenvalues
 * @param pvalues	p-values
 * @param twstat	statitic of tw
 * @param effectn	effective n (estimated)
 * @param percentage	percentage of variance
*/
void write_data_tracyWidom(char *output_file, int M, double *values,
                           double *pvalues, double *twstat, double *effectn,
                           double *percentage);

/**
 * calculate tracy widom  (twtable global variable  !!!)
 * 
 * @param values	eigenvalues
 * @param pvalues	p-values
 * @param twstat	statitic of tw
 * @param effectn	effective n (estimated)
 * @param M		number of non-zero eigenvalues
 */
void tw(double *values, double *pvalues, double *twstat, double *effectn,
        int M);

/**
 * calculate p-value from twstat (twtable global variable  !!!)
 *
 * @param lambda	statitic of tw
 *
 * @return p-values
 */
double twtest(double lambda);

/**
 * sort in decreasing order (rosetta stone)
 *
 * @param a	vector
 * @param n	size of a
 */
void insertion_sort(double *a, int n);

/**
 * remove zeros (< 1e-10) from a vector
 * 
 * @param values	vector
 * @param M		size
 */
void clean_zeros(double **values, int *M);

/**
 * clean from zeros and sort (decreasing order) a vector
 * 
 * @param values	vector
 * @param M		size
 */
void clean_sort(double **values, int *M);

#endif                          // TRACYWIDOM_H

/** @} */
