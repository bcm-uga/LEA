/** VCF conversion. Public native conversion signatures are unchanged. */
#ifndef VCF2GENO_H
#define VCF2GENO_H

#include <stddef.h>

/* Count samples in a validated header; release all resources before R errors. */
int nb_cols_vcf(char *file);

void vcf2geno(char *input_file, char *output_file, int *N, int *M,
              char *snp_file, char *removed_file, int *removed);

/* Internal wrapper interface: return 1 on success, 0 after cleanup on failure.
 * Diagnostic storage belongs to the caller; paths remain borrowed throughout. */
int vcf2geno_checked(char *input_file, char *output_file, int *N, int *M,
                     char *snp_file, char *removed_file, int *removed,
                     char *message, size_t message_size);

#endif
