/**
 * @addtogroup data_lfmm
 * @ingroup LFMM
 * @{
 * @file data_lfmm.h
 *
 * @brief set of functions to manage data
 */

#ifndef DATA_LFMM_H
#define DATA_LFMM_H

#include "../matrix/data.h"
#include "register_lfmm.h"
#include "lfmm_algo.h"

/**
 * simulate values from a normal distribution of
 * mean m_A and inverse covavariance inv_cov_A
 * into A (with alpha_R somewhere TODO) 
 * 
 * @param A     the output A realization (of size KxN)
 * @param m_A   the conditional A mean matrix (of size KxN)
 * @param inv_cov_A     the conditional A covariance matrix (of size KxK)
 * @param alpha_R       the inverse of the residual variance
 * @param K     the number of latent factors
 * @param N     the number of individuals
 * @param num_thrd      the number of processes used
 */
void rand_matrix(double *A, double *m_A, double *inv_cov_A, double alpha_R,
                 int K, int N, int num_thrd);

/** 
 * calculate inv_cov = inv(alphaR A %*% t(A) + diag(alpha))
 *
 * @param inv_cov  the output matrix (of size KxK)
 * @param alpha    the hyperparameters for A
 * @param alpha_R  the inverse of the residual variance
 * @param V     the A matrix (of size KxM)
 * @param K     a size 
 * @param M     a size
 * @param num_thrd      the number of processes used
 */
void create_inv_cov(double *inv_cov, double *alpha, double alpha_R,
                    double *A, int K, int M, int num_thrd);

/** calculate m = A * (R - B'*C))' if mode = 1 or m = A * (R - B'*C) if mode = °
 * @param A	matrix (of size KxM)
 * @param R	data matrix (of size NxM)
 * @param B	matrix (of size JxN)
 * @param C	matrix (of size JxM)
 * @param m	matrix (of size NxK)
 * @param M	size
 * @param N	size
 * @param J	size
 * @param K	size
 * @param num_thrd	number of processe used
 * @param mode	mode
 */
void create_m(double *A, float *R, double *B, double *C, double *m,
              int M, int N, int J, int K, int num_thrd, int mode);

/** calculates the quantiles of dist for prob
 *
 * @param dist  probability distribution
 * @param prob	probability for the quantiles
 * @param n	size of dist
 * @param p	size of prob
 * @param res	output quantiles
*/
void quantiles(double *dist, double *prob, int n, int p, double *res);

/**
 * calculate lambda criterion of a p.values distribution 
 * (Olivier's version with median over 41 values)
 *
 * @param p	p values set
 * @param n	size of p
 *
 * @return lambda
 */
double lambda(double *p, int n);

/**
 * convert pvalues into qvalues with Benjaminy-Hochberg approximation
 *
 * @param pvalues	pvalues table
 * @param qvalues	output qvalues table
 * @param n		size of the tables
 */
void pvalue_qvalue(double *pvalues, double *qvalues, int n);

/**
 * compute zscore in the dth column of beta
 * 
 * @param zscore	the output matrix of zscores (of size nxnD)
 * @param sum	the sum along the GS chain along for each individual (of size 2n)
 * @param sum2	the squared sum along the GS chain along for each individual (of size 2n)
 * @param n	the number of loci
 * @param cur	the length of the GS chain (without burnin)
 * @param D	number of lines of sum or sum2	
 */
void zscore_calc(double *zscore, double *sum, double *sum2, int n, int cur,
                 int D);

/**
 * add beta to the sum, called sum
 *
 * @param beta	the values to add (of size n)
 * @param sum	where to add beta (of size n)
 * @param n	the size of sum and beta
 */
void update_sum(double *beta, double *sum, int n);

/**
 * add squared beta to the sum, sum
 *
 * @param beta	the values to add (of size n)
 * @param sum	where to add beta (of size n)
 * @param n	the size of sum and beta
 */
void update_sum2(double *beta, double *sum2, int n);

/**
 * update deviance parameter
 *
 * @param deviance	the parameter to update
 * @param cur	the current length of the GS chain
 * @param var	the variance
 * @param thrd_m	the current mean
 */
void update_deviance(double *deviance, int cur, double var, double thrd_m);

/**
 * update mean parameter
 *
 * @param mean	the parameter to update (of size n)
 * @param beta	the beta parameter to update the mean (of size n)
 * @param n	the length of mean and beta
 * @param cur	the current length of the GS chain
 */
void update_mean(double *beta, double *mean, int n, int cur);

/**
 * from column d of C (of size NxnD), create Cpp (of size Nx2) with 1 in first column and
 * the column d of C in second column
 *
 * @param C	the covariable parameter (of size NxnD)
 * @param N	the number of individuals
 * @param nD	the number of covariables
 * @param Cpp	the output parameter (of size Nx2)
 * @param d	the column of C to copy (from 0)
 * @param all	if all, copy all columns of C into Cpp
 */
void modify_C(double *C, int N, int nD, double *Cpp, int d, int all);

/**
 * write deviance and DIC in file_data file
 *
 * @param file_data	the file where to register
 * @param deviance	the deviance parameter
 * @param DIC	the DIC parameter
 */
void write_DIC(char *file_data, double deviance, double DIC);

/**
 * compute and write absolute zscore, -log10(pvalue) and pvalue into file_data with %G 
 * separated by spaces
 *
 * @param file_data	the file where to write
 * @param M		the number of snps (of lines of dat) 
 * @param zscore	zscore matrix
 * @param D		D number of environmental variables
 * @param all		if all variables together
 * @param nd		number of the current env variable
 * @param K		number of latent factors
 * @param N		number of individuals
 * @param dev		deviance
 * @param DIC		DIC
 */
void write_zscore_double(char *output_file, int M, double *zscore, int D,
                         int all, int nd, int K, int N, double dev, double DIC);

/**
 * compute the current residual variance
 *
 * @param param		parameter structure
 * @param GS_param	GS parameter structure
 */
double var_data(LFMM_param param, LFMM_GS_param GS_param);

/**
 * compute the current residual variance and input missing data
 *
 * @param R	the data matrix (of size NxM)
 * @param I     the missing data matrix
 * @param U	the U matrix (of size KxN)
 * @param V	the V matrix (of size KxM)
 * @param C	the C matrix (of size NxD)
 * @param beta	the beta parameter (of size DxM)
 * @param N	the number of individuals
 * @param M	the number of loci
 * @param D	the number of covariables
 * @param K	the number of latent factors
 * @param thrd_m2	the output value of the mean
 * @param num_thrd	the number of thread used
 */
//double var_data_inputation(float *R, int *I, double *U, double *V, double *C, 
//      double *beta, int N, int M, int K, int D, double *thrd_m2, int num_thrd);

/**
 * input missing values with U'*V+C*beta
 * (not used)
 *
 * @param R     the data matrix (of size NxM)
 * @param U     the U matrix (of size DxN)
 * @param V     the V matrix (of size DxM)
 * @param C     the C matrix (of size NxK)
 * @param beta  the beta parameter (of size KxM)
 * @param I     the missing data matrix
 * @param N     the number of individuals
 * @param M     the number of loci
 * @param D     the number of covariables
 * @param K     the number of latent factors
 */
void inputation_lfmm(float *R, double *U, double *V, double *C, double *beta,
                     int *I, int N, int M, int K, int D);

/**
 * input missing values with empirical frequencies
 * (used for k=0)
 *
 * @param R     the data matrix (of size NxM)
 * @param I     the missing data matrix
 * @param N     the number of individuals
 * @param M     the number of loci
 */
void inputation_freq(float *R, int *I, int N, int M);

#endif                          // DATA_LFMM_H

/** @} */
