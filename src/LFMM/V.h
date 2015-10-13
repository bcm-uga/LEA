/**
 * @addtogroup V
 * @ingroup LFMM
 * @{
 * @file V.h
 *
 * @brief functions to update V parameters in LFMM model
 */

#ifndef V_H
#define V_H
#include "data_lfmm.h"

/**
 * update V matrix
 *
 * @param param         parameter structure
 * @param GS_param      GS parameter structure
 */
void update_V(LFMM_param param, LFMM_GS_param GS_param);

#endif                          // V_H

/** @} */
