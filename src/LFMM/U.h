/**
 * @addtogroup U
 * @ingroup LFMM
 * @{
 * @file U.h
 *
 * @brief functions to update U parameters in LFMM model
 */

#ifndef U_H
#define U_H
#include "data_lfmm.h"

/**
 * update U matrix
 *
 * @param param         parameter structure
 * @param GS_param      GS parameter structure
 */
void update_U(LFMM_param param, LFMM_GS_param GS_param);

/**
 * update alpha_U hyperparameter
 *
 * @param param         parameter structure
 */
void update_alpha_U(LFMM_param param);

#endif                          // U_H

/** @} */
