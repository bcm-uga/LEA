/**
 * @addtogroup slice_lfmm
 * @ingroup LFMM
 * @{
 * @file slice_lfmm.h
 *
 * @brief multithreaded part of the functions to compute 
 *	  new values for U, V or beta from the conditional 
 *	  distribution. 
 */

#ifndef SLICE_LFMM_H
#define SLICE_LFMM_H

#ifndef WIN32

/**
 * compute a slice of the lines of the conditional mean 
 *
 * @param G     a specific structure for multi-threading
 */
void slice_m(void *G);

/**
 * compute a slice of the columns of the new values 
 *
 * @param G     a specific structure for multi-threading
 */
void slice_rand(void *G);

/**
 * compute a slice of the lines of the conditional covariance
 *
 * @param G     a specific structure for multi-threading
 */
void slice_inv_cov(void *G);

#endif

#endif                          // SLICE_LFMM_H

/** @} */
