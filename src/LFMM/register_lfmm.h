/**
 * @addtogroup register_lfmm
 * @ingroup LFMM
 * @{
 * @file register_lfmm.h
 *
 * @brief functions to analyse the command line, init and and free the parameters. 
 */

#include "LFMM.h"

#ifndef REGISTER_LFMM_H
#define REGISTER_LFMM_H

/**
 * initialize the parameter with the default values for LFMM 
 *
 * @param param	parameter structure
 */
void init_param_lfmm(LFMM_param param);

/**
 * analyse command line set of parameters and set the parameters
 * 
 * @param argc	the number of arguments
 * @param argv	the set of arguments
 * @param param	parameter structure
 */
void analyse_param_lfmm(int argc, char *argv[], LFMM_param param);

/** 
 * free lfmm_param struct allocated memory
 *
 * @param param	parameter structure
 */
void free_param_lfmm(LFMM_param param);

#endif                          // REGISTER_LFMM_H

/** @} */
