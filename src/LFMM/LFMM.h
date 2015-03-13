/**
 * @addtogroup LFMM
 * @ingroup LFMM
 * @{
 * @file LFMM.h
 *
 * @brief functions to estimate the parameters of LFMM 
 *
 * The article describing the model is:
 * 
 * E. Frichot, S.D. Schoville, G. Bouchard, O. François. Testing for associations between loci and environmental gradients using latent factor mixed models. Molecular Biology and Evolution (2013) 30 (7): 1687-1699.
 * 
 * A supplementary material of the article describes the Gibbs Sampler algorithm.
 */

#ifndef LFMM_H
#define LFMM_H

/**
 * pointer to lfmm_param struct
 */
typedef struct _lfmm_param *LFMM_param;

/**
 * @brief Structure containing all parameters for LFMM
 */
typedef struct _lfmm_param {
        // data size parameters
        int D;                  /**< @brief the number of the variables in the variable file */
        int K;                  /**< @brief the number of latent factors */
        int nd;                 /**< @brief nd-th variable to use. It is possible to run LFMM only with 
				* the nd-th variable of the variable file */

        // Gibbs Sampler parameters
        int Niter;              /**< @brief the total number of iterations (with burnin)*/
        int burn;               /**< @brief the number of burnin iterations*/
        int num_thrd;           /**< @brief the number of processes used */
        int init;               /**< @brief if true, random init. Otherwise, init with zeros */

        // model parameters
        double *alpha_beta;     /**< @brief the vector of hyperparameters for beta (of size D) */
        double alpha_R;         /**< @brief the inverse of the residual variance */
        double *alpha_U;        /**< @brief the hyperparameter for U */
        double *alpha_V;        /**< @brief the hyperparameter for V */
        double noise_epsilon;   /**< @brief prior for the different variances */
        double b_epsilon;       /**< @brief prior for the correlation coefficient variance */
        int mD;                 /**< @brief the number of the variables for the LFMM run (including intercept)*/

        // init parameters
        int *I;                 /**< @brief missing data matrix */
        int missing_data;       /**< @brief boolean: true if missing data */
        long long seed;         /**< @brief seed values */
        int all;                /**< @brief Boolean, if true, all variables are included at the same time 
				 * in a LFMM run. If false, each varriable is included in a separate run */

        // matrix parameters
        double *U;              /**< @brief the U matrix (of size Kxn)*/
        double *V;              /**< @brief the V matrix (of size KxL)*/
        float *dat;             /**< @brief the data matrix (of size nxL) */
        double *beta;           /**< @brief the beta matrix (of size mDxL)*/
        double *C;              /**< @brief The variable matrix (of size nxD), read from the variable file */
        double *mC;             /**< @brief The variable matrix for LFMM run (of size mDxn (warning: the size
				 * is different from C)) */
        double *zscore;         /**< @brief the output zscore matrix */

        // io parameters
        char output_file[512];  /**< @brief output file */
        char input_file[512];   /**< @brief input file */
        char cov_file[512];     /**< @brief variable file */
        int n;                  /**< @brief the number of individuals */
        int L;                  /**< @brief the number of loci */

        // output criterion parameters
        double dev;             /**< @brief deviance criterion */
        double DIC;             /**< @brief DIC criterion */
} lfmm_param;

/**
 * run LFMM
 *
 * @param param		parameter structure
 */
void LFMM(LFMM_param param);

#endif                          // LFMM_H

/** @} */
