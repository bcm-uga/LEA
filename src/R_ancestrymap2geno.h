/**
 * @file R_ancestrymap2geno.h
 *
 * @brief C wrapper for ancestrymap2geno command line program
 */

#ifndef R_ancestrymap2geno_H
#define R_ancestrymap2geno_H

void R_ancestrymap2geno(char **R_input_file, char **R_output_file, int *N,
                        int *M);

#endif                          // R_ancestrymap2geno_H
