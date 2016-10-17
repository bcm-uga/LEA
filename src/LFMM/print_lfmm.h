/**
 * @addtogroup print_lfmm
 * @ingroup LFMM
 * @{
 * @file print_lfmm.h
 *
 * @brief set of printing functions
 */

#ifndef PRINT_H
#define PRINT_H

#include "register_lfmm.h"

/**
 * print the complete licence
 */
void print_licence_lfmm();

/**
 * print the header for the licence
 */
void print_head_licence_lfmm();

/**
 * print my header
 */
void print_head_lfmm();

/**
 * print the help
 */
void print_help_lfmm();

/**
 * print summary of the parameters
 *
 * @param param	parameter structure
 */
void print_summary_lfmm(LFMM_param param);

/**
 * print command line options 
 *
 * @param argc	the number of options 
 * @param argv	the set of arguments
 */
void print_options_lfmm(int argc, char *argv[]);

/**
 * print percentage of variance 
 *
 * @param perc	percentage of variance
 * @param K	number of latent factors
 * @param D	number of variables
 */
void print_perc(double *perc, int K, int D);

#endif                          // PRINT_H

/** @} */
