/**
 * @addtogroup register_tracyWidom
 * @ingroup tracyWidom
 * @{
 * @file register_tracyWidom.h
 *
 * @brief read cmd line parameter for tracyWidom program
 */

#ifndef REGISTER_TRACYWIDOM_H
#define REGISTER_TRACYWIDOM_H

/**
 * analyse command line set of parameters and set the parameters
 * 
 * @param argc  	the number of arguments
 * @param argv  	the set of arguments
 * @param input		the input file
 * @param output	the output file	
 */
void analyse_param_tracyWidom(int argc, char *argv[], char *input,
                              char *output);

#endif                          // REGISTER_TRACYWIDOM_H

/** @} */
