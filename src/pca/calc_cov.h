/**
 * @addtogroup calc_cov
 * @ingroup pca
 * @{
 * @file calc_cov.h
 *
 * @brief functions to calculate the covariance matrix,
 *	x and the standard deviation. 
 */

#ifndef CALC_COV_H
#define CALC_COV_H

/**
 * calculate the covariance matrix of data into cov 
 *
 * @param data		data matrix (of size NxM)
 * @param N		number of individuals
 * @param M		number of loci
 * @param cov		output covariance matrix
 */
void calc_cov(double *data, int N, int M, double *cov);

/**
 * calculate sdev form eigenvalues, sdev = sqrt(val/N)
 *
 * @param val	eigenvalues vector
 * @param N	size of val
 */
void calc_sdev(double *val, int N);

/**
 * calculate x from vect and val (in fact sdev), vec <-x = vec * val /sqrt(N)
 *
 * @param vec	in, vec matrix. out, x (of size NxK)
 * @param val	sdev matrix
 * @param N	number of individuals
 * @param K	number of PCs
 */
void calc_x(double *vec, double *val, int N, int K);

#endif // CALC_COV_H

/** @} */
