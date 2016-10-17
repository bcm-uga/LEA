/**
 * @addtogroup lfmm_algo
 * @ingroup LFMM
 * @{
 * @file lfmm_algo.h
 *
 * @brief functions to compute lfmm parameters (zscores)
 * 
 */

#ifndef lfmm_ALGO_H
#define lfmm_ALGO_H

#include "register_lfmm.h"

/**
 * pointer to lfmm_GS_param struct
 */
typedef struct _lfmm_GS_param *LFMM_GS_param;

/** 
 * @brief structure containing paramaters for the LFMM Gibbs Sampler
 */
typedef struct _lfmm_GS_param {
        // beta
        double *m_beta;         /**< @brief conditional mean for beta (DxL)*/
        double *inv_cov_beta;   /**< @brief conditional inverse covariance matrix for beta */
        double *sum;            /**< @brief sum value for b*/
        double *sum2;           /**< @brief sum of square values of b*/
        // U
        double *m_U;            /**< @brief conditional mean for U (Kxn)*/
        double *mean_U;         /**< @brief mean value for U*/
        double *inv_cov_U;      /**< @brief conditional inverse covariance matrix for U */
        // V
        double *m_V;            /**< @brief conditional mean for V (KxL)*/
        double *mean_V;         /**< @brief mean value for V*/
        double *inv_cov_V;      /**< @brief conditional inverse covariance matrix for V */

        // ?
        double thrd_m2;         /**< @brief It contains ||G-UV^T -XB^T||_F^2. It is updated at each
				 * iteration in var_data through update_alpha_R.
				 * It is used to calculate the deviance and the DIC criterion */

} lfmm_GS_param;

/**
 * calculate LFMM parameters
 *
 * @param param	parameter structure
 */
void lfmm_emcmc(LFMM_param param);

/**
 * update different sums
 *
 * @param param		parameters structure 
 * @param GS_param	GS parameters structure 
 */
void update_sums(LFMM_param param, LFMM_GS_param GS_param);

/**
 * update dp and deviance
 *
 * @param param		parameters structure 
 * @param GS_param	GS parameters structure 
 * @param deviance	deviance
 * @param dp		temporary calc
 */
void calc_dp_deviance(LFMM_param param, LFMM_GS_param GS_param,
                      double *deviance, double *dp);

/**
 * allocate all temporary memory
 *
 * @param GS_param	GS parameters structure 
 * @param N     the number of individuals
 * @param M     the number of loci
 * @param K     the number of latent factors
 * @param D     the number of covariables
 */
void allocate_all(LFMM_GS_param GS_param, int N, int M, int K, int D);

/**
 * free all temporary memory
 *
 * @param GS_param	GS parameters structure 
 */
void free_all(LFMM_GS_param GS_param);

/** 
 * update alpha_R parameter 
 *
 * @param param		parameters structure 
 * @param GS_param	GS parameters structure 
 */
double update_alpha_R(LFMM_param param, LFMM_GS_param GS_param);

#endif                          // lfmm_algo_H

/** @} */
