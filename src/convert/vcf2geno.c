/*
 *     vcf2geno, file: vcf2geno.c
 *     Copyright (C) 2013 Eric Frichot
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
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
#include "../io/io_tools.h"
#include "../io/io_error.h"
#include "../io/read.h"
#include "geno.h"
#include "vcf2geno.h"

// read_lines

int nb_cols_vcf(char *file)
{
        char tok;
        FILE *File = fopen_read(file);
        int max = 1, cols = 0;
        int N = 0;
        char *szbuff;
        char *token;

        tok = fgetc(File);
        while (!feof(File) && tok == '#') {
                // count number of elements
                cols = 1;
                while (tok != '\n' && tok != EOF) {
                        tok = fgetc(File);
                        cols++;
                }
                // save the max
                if (cols > max)
                        max = cols;
                if (!feof(File))
                        tok = fgetc(File);
        }
        fclose(File);

        File = fopen_read(file);
        szbuff = (char *)Calloc(2 * max * sizeof(char), char);
        token = fgets(szbuff, 2 * max, File);

        // skip comment lines
        while (szbuff[1] == '#')
                token = fgets(szbuff, 2 * max, File);

        // count the number of elements
        token = strtok(szbuff, SEP);
        if (!strcmp(token, "#CHROM")) {
                while (token) {
                        N++;
                        token = strtok(NULL, SEP);
                }
        } else {
                Rprintf
                    ("Error: It seems that the line with individual informations is missing.\n");
                error(NULL);
        }

        // close and free
        fclose(File);
        Free(szbuff);

        // remove the 9 information elements
        return N - 9;
}

// vcf2geno

void vcf2geno(char *input_file, char *output_file, int *N, int *M,
              char *snp_file, char *removed_file, int *removed)
{
        // temporary variables
        int i, j, jp, jt;
        char **infos = (char **)Calloc(9 *  sizeof(char *), char *);

        // file management
        FILE *input_File = NULL;
        FILE *snp_File = NULL;
        FILE *removed_File = NULL;
        FILE *output_File = NULL;
        int max_char_per_line;
        char *szbuff, *token = NULL;

        // register output
        int *allele;

        // read comments lines and first lines
        *N = nb_cols_vcf(input_file);
        if ((*N) <= 0) {
                Rprintf("Error: It seems that %s (vcf file) contains no genotype"
                       " information.\n", input_file);
                error(NULL);
        }

        // init tmp mem
        for (i = 0; i < 9; i++)
                infos[i] = (char *)Calloc(512 * sizeof(char), char);

        max_char_per_line = 50 * (*N) + 20;
        szbuff = (char *)Calloc(max_char_per_line * sizeof(char), char);
        allele = (int *)Calloc(*N * sizeof(int), int);

        // open files
        input_File = fopen_read(input_file);
        snp_File = fopen_write(snp_file);
        removed_File = fopen_write(removed_file);
        output_File = fopen_write(output_file);

        j = 0;
        jp = 0;
        jt = 0;

        // for each line
        while (fgets(szbuff, max_char_per_line, input_File)) {
                // if not a comment line
                if (szbuff[0] != '#') {

                        // read CNV informations
                        read_cnv_info(token, infos, szbuff, jt + 1);

                        // if not a SNP.
                        if (strlen(infos[3]) > 1 || strlen(infos[4]) > 1) {
                                // print it in removed lines info file as REMOVED
                                write_snp_info(removed_File, infos, 1);
                                jp++;

                                // if it is a SNP
                        } else {
                                // print it in SNP info file
                                write_snp_info(snp_File, infos, 0);

                                // read line vcf
                                fill_line_vcf(token, allele, j, *N,
                                              input_file, input_File);

                                // print geno line
                                write_geno_line(output_File, allele, *N);
                                j++;
                        }
                }
                // next line
                jt++;
        }
        *removed = jp;
        *M = j;

        // close files
        fclose(input_File);
        fclose(snp_File);
        fclose(removed_File);
        fclose(output_File);

        // free memory
        for (i = 0; i < 9; i++)
                Free(infos[i]);
        Free(infos);
        Free(szbuff);
        Free(allele);
}

// fill_line_vcf

void fill_line_vcf(char *token, int *allele, int j, int N, char *input_file,
                   FILE * input_File)
{

        int i = 0;

        // for each otoken
        token = strtok(NULL, SEP);
        while (token) {
                // if length 1
                if (strlen(token) == 1) {
                        // and just a dot
                        if (token[0] == '.')
                                allele[i] = 9;
                        else {
                                Rprintf
                                    ("Error: SNP %d, individual %d, not 0/1.\n",
                                     j + 1, i + 1);
                                error(NULL);
                        }
                        // if just a dot and a EOL
                } else if (strlen(token) == 2 && token[0] == '.'
                           && token[1] == '\n') {
                        allele[i] = 9;
                        // if not [01]/[01]
                } else if (token[1] == ':') {
                        Rprintf("Error: SNP %d, individual %d, not a SNP.\n",
                               j + 1, i + 1);
                        error(NULL);
                } else {
                        // if error
                        if ((token[0] != '0' && token[0] != '1'
                             && token[0] != '.')
                            || (token[2] != '0' && token[2] != '1'
                                && token[2] != '.')) {
                                Rprintf
                                    ("Error: SNP %d, individual %d, not a 0/1: '%c','%c'.\n",
                                     j + 1, i + 1, token[0], token[2]);
                                error(NULL);
                        }
                        // if missing data
                        if (token[0] == '.' || token[2] == '.')
                                allele[i] = 9;
                        // else everything is fine
                        else {
                                allele[i] =
                                    (int)token[0] - (int)'0' + (int)token[2] -
                                    (int)'0';
                        }
                }
                // next token
                token = strtok(NULL, SEP);
                i++;

        }

        // if not the correct number of columns
        test_column(input_file, input_File, i, j + 1, N, token);
}

// read_cnv_info

void read_cnv_info(char *token, char **infos, char *szbuff, int j)
{
        int ip;
        token = strtok(szbuff, SEP);
        if (!token) {
                Rprintf("Error while reading SNPs informations at line %d.\n",
                       j);
                error(NULL);
        }

        strcpy(infos[0], token);
        for (ip = 1; ip < 9; ip++) {
                token = strtok(NULL, SEP);
                if (!token) {
                        Rprintf
                            ("Error while reading SNPs informations at line %d.\n",
                             j);
                        error(NULL);
                }
                strcpy(infos[ip], token);
        }
}

// write_snp_info

void write_snp_info(FILE * output_File, char **infos, int removed)
{
        int ip;
        for (ip = 0; ip < 9; ip++)
                fprintf(output_File, "%s ", infos[ip]);

        if (removed)
                fprintf(output_File, "REMOVED\n");
        else
                fprintf(output_File, "\n");
}
