/**
 * @addtogroup beta_k0
 * @ingroup LFMM
 * @{
 * @file beta_k0.h
 *
 * @brief functions to calc beta parameters in LM model.
 * 	The code is based on a previous matlab version. 
 *	Comments in matlab code format can be found in the code.
 */

#ifndef BETA_K0_H
#define BETA_K0_H
#include "data_lfmm.h"

/**
 * calc beta matrix
 *
 * @param C	the covariable matrix (of size NxK)
 * @param data	the data matrix (of size NxM)
 * @param beta	the output beta realization (of size KxM)
 * @param CCt	the constant covariance matrix for C
 * @param var_beta	the conditional beta mean matrix (of size KxM)
 * @param M	the number of loci
 * @param N	the number of individuals
 * @param D	the number of covariables
 * @param var_res	variance of the residual	
 */
void calc_beta_k0(double *C, float *R, double *beta,
                  double *CCt, double *var_beta, int M, int N,
                  int D, double *var_res);

/**
 * calculate zscores 
 * 
 * @param zscore	output zscore
 * @param beta		esperance of beta
 * @param var_beta	variance of beta
 * @param D		number of lines of each
 * @param M		number of columns of each
 */
void zscore_calc_k0(double *zscore, double *beta, double *var_beta, int D,
                    int M);

/**
 * compute the C'*C in tmp
 *
 * @param cov   the output matrix (of size DxD)
 * @param C     the covariable matrix (of size NxD)
 * @param D     the number of covariables
 * @param N     the number of individuals
 */
void create_CCt(double *cov, double *C, int D, int N);

#endif                          // BETA_K0_H

/** @} */
