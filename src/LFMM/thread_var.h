/**
 * @addtogroup thread_var
 * @ingroup LFMM
 * @{
 * @file thread_var.h
 *
 * @brief functions to compute new values for the variance
 *        (possibly multithreaded) 
 */

#ifndef THREAD_VAR_H
#define THREAD_VAR_H

#ifndef WIN32

#include "register_lfmm.h"
#include "lfmm_algo.h"

/**
 * structure to manage multithreading to calculate the variance and the mean
 */
typedef struct _multithreading_lfmm_var *Multithreading_lfmm_var;

/**
 * @brief Structure containing generic parameters for multithreading variance 
 *	calculation in LFMM functions. 
 */
typedef struct _multithreading_lfmm_var {
        float *R;
        double *U;
        double *V;
        double *C;
        double *beta;
        int D;
        int N;
        int M;
        int K;
        double mean;
        double res;
        double res2;
        int slice;
        int num_thrd;
} multithreading_lfmm_var;

/**
 * general multithreading function manager. Some parameters can be NULL
 *
 * @param param		parameter structure
 * @param GS_param	GS parameter structure
 * @param fct		generic pointer function
 * @param res	the first res
 * @param res2	the second res
 */
void thrd_var(LFMM_param param, LFMM_GS_param GS_param,
              void (*fct) (), double *res, double *res2);

/**
 * compute a slice of the mean
 *
 * @param G     a specific structure for multi-threading
 */
void slice_mean(void *G);

/**
 * compute a slice of the var
 *
 * @param G     a specific structure for multi-threading
 */
void slice_var(void *G);

#endif                          // THREAD_VAR_H

/** @} */

#endif
