/**
 * @addtogroup print_tracyWidom
 * @ingroup tracyWidom
 * @{
 * @file print_tracyWidom.h
 *
 * @brief functions to print help and a summary of the parameters
 */

#ifndef PRINT_TRACYWIDOM_H
#define PRINT_TRACYWIDOM_H

/**
 * print help
 */
void print_help_tracyWidom();

/**
 * print summary of the parameters
 *
 * @param M     		the number of eigenvalues
 * @param input         	input eigenvalues file
 * @param output 	  	output file
 */
void print_summary_tracyWidom(int M, char *input, char *output);

#endif                          // PRINT_TRACYWIDOM_H

/** @} */
