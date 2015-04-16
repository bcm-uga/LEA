/*
   vcf2geno, file: main.c
   Copyright (C) 2013 François Mathieu, Eric Frichot

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
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
#include <time.h>
#include <math.h>
#include "convert/vcf2geno.h"
#include "io/io_tools.h"
#include "convert/register_convert.h"

#include "R_vcf2geno.h"

void R_vcf2geno(char **R_input_file, char **R_output_file, int *N, int *M)
{
        int removed;            // number of removed CNVs
        char snp_bp_file[512];  // snp file
        char removed_bp_file[512];      // removed snp file
        char *tmp;

        tmp = remove_ext(*R_output_file, '.', '/');
        strcpy(snp_bp_file, tmp);
        strcat(snp_bp_file, ".vcfsnp");

        strcpy(removed_bp_file, tmp);
        strcat(removed_bp_file, ".removed");

        vcf2geno(*R_input_file, *R_output_file, N, M, snp_bp_file,
                 removed_bp_file, &removed);

        print_convert(*N, *M);

        Rprintf("For SNP info, please check %s.\n\n", snp_bp_file);

        Rprintf("%d line(s) were removed because these are not SNPs.\n",
                removed);
        Rprintf("Please, check %s file, for more informations.\n\n",
                removed_bp_file);

        free(tmp);

}

        /*
           // parameter initialization
           char input_file[512];                // input file "without" missing data
           char output_file[512];                // output file with missing data
           double e = 0.05;                // output percentage of missing data
           long long seed = -1;
           int m = 0;

           print_head_snmf();

           if (R_input_file)
           strcpy(input_file, *R_input_file);
           else 
           print_error_cds("option","-g vcftype_file");

           createDataSet(input_file, m, (long long) seed, percentage, output_file);
         */
