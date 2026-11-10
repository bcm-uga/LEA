/*
 *     convert, file: ancestrymap.c
 *     Copyright (C) 2013 Eric Frichot
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <R.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ancestrymap.h"
#include "geno.h"
#include "../io/io_tools.h"
#include "../io/io_error.h"
#include "../io/io_data_int.h"
#include "../io/read.h"
#include "register_convert.h"

/* Validate before either SNP-name copy. Lengths are bytes, excluding NUL. */
static void check_snp_name(const char *token, const char *input, int i, int j)
{
        if (!token)
                error("Invalid ANCESTRYMAP record in %s at individual %d, SNP %d: "
                      "missing SNP name.", input, i, j);
        if (strlen(token) >= ANCESTRYMAP_SNP_NAME_CAPACITY)
                error("SNP name too long in %s at individual %d, SNP %d: "
                      "maximum is %d bytes.", input, i, j,
                      ANCESTRYMAP_SNP_NAME_CAPACITY - 1);
}

// ancestrymap2geno

void ancestrymap2geno(char *input_file, char *output_file, int *N, int *M)
{
        int *data;
        int nb;
        double tmp;

        // number of lines and columns
        *N = nb_ind_ancestrymap(input_file);
        nb = nb_lines(input_file, 1000);
        tmp = (double)nb / (double)(*N);
        if (tmp != floor(tmp)) {
                error("Invalid ANCESTRYMAP file %s: incorrect number of lines "
                      "for %d individuals.", input_file, *N);
        }
        *M = (int)tmp;

        // allocate_memory
        data = (int *) calloc((*N) * (*M) ,  sizeof(int));

        // read in ancestrymap format
        read_ancestrymap(input_file, *N, *M, data);

        // write in geno format
        write_geno(output_file, *N, *M, data);

        // free memory
        free(data);
}

// ancestsrymap2lfmm

void ancestrymap2lfmm(char *input_file, char *output_file, int *N, int *M)
{
        int *data;
        int nb;
        double tmp;

        // set the number of individuals (N) and the number of loci (M)
        *N = nb_ind_ancestrymap(input_file);
        nb = nb_lines(input_file, 1000);
        tmp = (double)nb / (double)(*N);
        if (tmp != floor(tmp)) {
                error("Invalid ANCESTRYMAP file %s: incorrect number of lines "
                      "for %d individuals.", input_file, *N);
        }
        *M = (int)tmp;

        // allocate memory
        data = (int *) calloc((*N) * (*M) , sizeof(int));

        // read in ancestrymap format
        read_ancestrymap(input_file, *N, *M, data);

        // write in lfmm format
        write_data_int(output_file, *N, *M, data);

        // free memory
        free(data);
}

// read_ancestrymap

void read_ancestrymap(char *input_file, int N, int M, int *data)
{
        FILE *File = NULL;
        int i, j, allele;
        char tmp[ANCESTRYMAP_SNP_NAME_CAPACITY] = "";
        char ref[ANCESTRYMAP_SNP_NAME_CAPACITY] = "";
        int max_char_per_line = 1000;
        char szbuff[max_char_per_line];
        int warning = 0;

        // open input file
        File = fopen_read(input_file);

        j = 0;
        i = 0;
        // while not the end of the file
        //      or not too many SNPs (j < M-1)
        //      or last SNP but not too many ind (j == (M-1) && i< N))
        while (fgets(szbuff, max_char_per_line, File)
               && (j < M - 1 || (j == (M - 1) && i < N))) {

                // read line
                read_line_ancestrymap(szbuff, &allele, tmp, i + 1, j + 1,
                                      input_file, &warning);

                // first SNP, save SNP name
                if (j == 0 && i == 0)
                        strcpy(ref, tmp);

                // new SNP
                if (strcmp(ref, tmp) != 0) {
                        // test if the number individual for the line is ok
                        if (i != N) {
                                fclose(File);
                                error("Invalid ANCESTRYMAP file %s at SNP %d: "
                                      "%d individuals instead of %d.",
                                      input_file, j + 1, i, N);
                        }
                        // new line
                        i = 0;
                        j++;
                        // save new SNP name
                        strcpy(ref, tmp);
                }
                // Reject inconsistent groups before indexing the data matrix.
                if (i >= N || j >= M) {
                        fclose(File);
                        error("Invalid ANCESTRYMAP group sizes in %s at "
                              "individual %d, SNP %d.", input_file, i + 1, j + 1);
                }
                // write genotype
                data[i * M + j] = allele;

                i++;
        }
        if (ferror(File)) {
                fclose(File);
                error("Error reading ANCESTRYMAP file %s.", input_file);
        }
        if (i != N || j + 1 != M || !feof(File)) {
                fclose(File);
                error("Invalid ANCESTRYMAP file %s: inconsistent number of "
                      "individuals or loci.", input_file);
        }

        fclose(File);
}

// read_line

void read_line_ancestrymap(char *szbuff, int *allele, char *name, int i, int j,
                           char *input, int *warning)
{
        char *token;

        // read and validate the SNP name before copying into the caller's buffer
        token = strtok(szbuff, SEP);
        check_snp_name(token, input, i, j);
        strcpy(name, token);

        // The individual identifier and genotype columns must both exist.
        token = strtok(NULL, SEP);
        if (!token)
                error("Invalid ANCESTRYMAP record in %s at individual %d, SNP %d: "
                      "missing individual identifier and genotype.", input, i, j);
        token = strtok(NULL, SEP);
        if (!token)
                error("Invalid ANCESTRYMAP record in %s at individual %d, SNP %d: "
                      "missing genotype.", input, i, j);
        *allele = (int)atoi(token);

        // if not 0, 1, 2, or 9
        if (!(*warning) && *allele != 9 && *allele != 1 && *allele != 2
            && *allele != 0) {
                Rprintf("Warning: some genotypes are not 0, 1, 2 or 9 in %s.\n",
                       input);
                Rprintf("\t First warning at individual %d, column %d.\n\n", j,
                       i);
                *warning = 1;
        }
}

// nb_ind_ancestrymap

int nb_ind_ancestrymap(char *input_file)
{
        FILE *File = NULL;
        int max_char_per_line = 1000;
        int nb = 0;
        char szbuff[max_char_per_line];
        char *token;
        char tmp[ANCESTRYMAP_SNP_NAME_CAPACITY] = "";

        File = fopen_read(input_file);

        // Count the first locus, using only successfully read records.
        while (fgets(szbuff, max_char_per_line, File)) {
                token = strtok(szbuff, SEP);
                check_snp_name(token, input_file, nb + 1, 1);
                if (nb == 0)
                        strcpy(tmp, token);
                else if (strcmp(tmp, token) != 0)
                        break;
                nb++;
        }
        if (ferror(File)) {
                fclose(File);
                error("Error reading ANCESTRYMAP file %s.", input_file);
        }
        fclose(File);
        if (nb == 0)
                error("Empty ANCESTRYMAP file %s: no records.", input_file);
        return nb;
}
