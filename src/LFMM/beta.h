/**
 * @addtogroup beta
 * @ingroup LFMM
 * @{
 * @file beta.h
 *
 * @brief functions to update beta parameters in LFMM model.
 * 	The code is based on a previous matlab version. 
 *	Comments in matlab code format can be found in the code.
 */

#ifndef BETA_H
#define BETA_H
#include "data_lfmm.h"
#include "lfmm_algo.h"
#include "register_lfmm.h"

/**
 * sample new values for the beta matrix given its conditional distribution
 *
 * @param param 	parameter structure
 * @param GS_param 	GS parameter structure 
 */
void update_beta(LFMM_param param, LFMM_GS_param GS_param);

/**
 * sample new values for the alpha_beta hyperparameters given its conditional 
 * distribution
 *
 * @param param 	parameter structure
 */
void update_alpha_beta(LFMM_param param);

#endif                          // BETA_H

/** @} */
