/**
 * @file R_vcf2geno.h
 *
 * @brief C wrapper for vcf2geno command line program
 */

#ifndef R_vcf2geno_H
#define R_vcf2geno_H

void R_vcf2geno(char **R_input_file, char **R_output_file, int *N, int *M);

#endif                          // R_vcf2geno_H
