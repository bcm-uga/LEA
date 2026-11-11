/*
   vcf2geno, file: main.c
   Copyright (C) 2013 François Mathieu, Eric Frichot

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "convert/vcf2geno.h"
#include "convert/register_convert.h"
#include "R_vcf2geno.h"

void R_vcf2geno(char **R_input_file, char **R_output_file, int *N, int *M)
{
        int removed, success;
        const char *output, *dot, *slash;
        size_t length;
        char *snp_file, *removed_file;
        char message[1024];

        if (!R_input_file || !*R_input_file || !R_output_file || !*R_output_file)
                error("VCF conversion requires input and output filenames.");
        output = *R_output_file;
        dot = strrchr(output, '.');
        slash = strrchr(output, '/');
        /* Match remove_ext(output, '.', '/') without an extra owned copy. */
        length = dot && (!slash || slash < dot) ? (size_t)(dot - output) : strlen(output);
        if (length > SIZE_MAX - sizeof(".removed"))
                error("VCF output path exceeds allocation limits.");
        snp_file = malloc(length + sizeof(".vcfsnp"));
        removed_file = malloc(length + sizeof(".removed"));
        if (!snp_file || !removed_file) {
                free(snp_file);
                free(removed_file);
                error("Unable to allocate VCF output paths.");
        }
        memcpy(snp_file, output, length);
        memcpy(snp_file + length, ".vcfsnp", sizeof(".vcfsnp"));
        memcpy(removed_file, output, length);
        memcpy(removed_file + length, ".removed", sizeof(".removed"));

        success = vcf2geno_checked(*R_input_file, *R_output_file, N, M,
                                   snp_file, removed_file, &removed,
                                   message, sizeof(message));
        if (success) {
                print_convert(*N, *M);
                Rprintf("For SNP info, please check %s.\n\n", snp_file);
                Rprintf("%d line(s) were removed because these are not SNPs.\n", removed);
                Rprintf("Please, check %s file, for more informations.\n\n", removed_file);
        }
        free(snp_file);
        free(removed_file);
        if (!success)
                error("%s", message);
}
