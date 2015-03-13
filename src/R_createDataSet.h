/**
 * @file R_createDataSet.h
 *
 * @brief C wrapper for createDataSet command line program
 */

#ifndef R_createDataSet_H
#define R_createDataSet_H

void R_createDataSet(char **R_input_file, int *R_seed, double *R_percentage,
                     char **R_output_file);

#endif                          // R_createDataSet_H
