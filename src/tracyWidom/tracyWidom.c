/*
 * tracyWidom, file: tracyWidom.c Copyright (C) 2013 François Mathieu, Eric
 * Frichot
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with 
 * this program.  If not, see <http://www.gnu.org/licenses/>. 
 */
#include <R.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../io/io_tools.h"
#include "../io/io_data_double.h"
#include "tracyWidom.h"
#include "twtable.h"
#include "print_tracyWidom.h"

// tracyWidom

void tracyWidom(char *input_file, char *output_file)
{
        double *values, *pvalues, *twstat, *effectn, *percentage;
        int i, N, M;
        double sum;

        // number of lines and columns
        M = nb_cols_lfmm(input_file);
        N = nb_lines(input_file, M);

        // correct K
        if (M != 1) {
                Rprintf("Tracy-Widom: Error %s has more than one column\n",
                       input_file);
        }
        // print command line summary
        print_summary_tracyWidom(N, input_file, output_file);

        // allocate memory 
        values = (double *)Calloc(N *  sizeof(double), double);

        // read input_file
        read_data_double(input_file, N, 1, values);

        // sort value and remove 0 
        clean_sort(&values, &N);

        // allocate memory
        pvalues = (double *)Calloc(N *  sizeof(double), double);
        twstat = (double *)Calloc(N *  sizeof(double), double);
        effectn = (double *)Calloc(N *  sizeof(double), double);
        percentage = (double *)Calloc(N *  sizeof(double), double);

        // calculate tracy-widom values
        tw(values, pvalues, twstat, effectn, N);

        sum = 0.0;
        for (i = 0; i < N; i++)
                sum += values[i];
        for (i = 0; i < N; i++)
                percentage[i] = values[i] / sum;

        // write output
        write_data_tracyWidom(output_file, N, values, pvalues, twstat, effectn,
                              percentage);

        // free memory
        Free(values);
        Free(pvalues);
        Free(twstat);
        Free(effectn);
        Free(percentage);
}

// write_data_tracyWidom

void write_data_tracyWidom(char *output_file, int M, double *values,
                           double *pvalues, double *twstat, double *effectn,
                           double *percentage)
{
        FILE *file = NULL;
        int i;

        // open file
        file = fopen_write(output_file);

        // header
        fprintf(file,
                "N\teigenvalues\ttwstats\t\tpvalues\teffectn\tpercentage\n");
        for (i = 0; i < M; i++) {
                fprintf(file,
                        "%d\t%3.4G\t\t%3.4G\t\t%3.4G\t\t%3.8G\t%3.4G\n",
                        i + 1, values[i], twstat[i], pvalues[i], effectn[i],
                        percentage[i]);
        }

        fclose(file);
}

// tw

void tw(double *values, double *pvalues, double *twstat, double *effectn, int M)
{
        int Mp, i;
        double nhat, lambda, mu, sigma, s, s2;

        // s & s2
        s = 0;
        s2 = 0;
        for (i = 0; i < M; i++) {
                s += values[i];
                s2 += values[i] * values[i];
        }

        for (i = 0; i < M; i++) {
                // Mp
                Mp = M - i;
                // nhat
                nhat = (Mp + 2) * s * s / (Mp * s2 - s * s);
                // lambda
                lambda = values[i] * Mp / s;
                // mu
                mu = (sqrt(nhat - 1) + sqrt(Mp)) * (sqrt(nhat - 1) +
                                                    sqrt(Mp)) / nhat;
                // sigma
                sigma = (sqrt(nhat - 1) + sqrt(Mp)) / nhat
                    * pow(1.0 / sqrt(nhat - 1) + 1.0 / sqrt(Mp),
                          (double)(1.0 / 3.0));
                // lambda
                lambda = (lambda - mu) / sigma;
                // test tracyWidom
                pvalues[i] = twtest(lambda);
                twstat[i] = lambda;
                effectn[i] = nhat;
                // udpate s & s2
                s -= values[i];
                s2 -= values[i] * values[i];
        }
}

// twtest

double twtest(double lambda)
{
        int i = 0;

        // find the correct line
        while (i < 161 && lambda >= twtable[i * 3])
                i++;

        if (i == 161)
                return (twtable[(i - 1) * 3 + 1]);
        else if (i == 0)
                return (twtable[i * 3 + 1]);
        else {
                return (twtable[(i - 1) * 3 + 1] +
                        (twtable[i * 3 + 1] - twtable[(i - 1) * 3 + 1])
                        * (lambda - twtable[(i - 1) * 3])
                        / (twtable[i * 3] - twtable[(i - 1) * 3]));
        }

}

// insertion_sort (rosetta code)

void insertion_sort(double *a, int n)
{
        int i, j;
        double value;
        for (i = 1; i < n; i++) {
                value = a[i];
                for (j = i; j > 0 && value > a[j - 1]; j--) {
                        a[j] = a[j - 1];
                }
                a[j] = value;
        }
}

// clean_zeros

void clean_zeros(double **values, int *M)
{
        double *new, *tmp;
        int i;

        // find 0
        i = *M - 1;
        while (fabs((*values)[i]) < 1e-10)
                i--;

        // modify 
        i++;
        if (i < *M) {
                *M = i;
                new = (double *)Calloc(*M *  sizeof(double), double);
                for (i = 0; i < *M; i++)
                        new[i] = (*values)[i];
                tmp = *values;
                *values = new;
                Free(tmp);
        }
}

// clean_sort

void clean_sort(double **values, int *M)
{
        // sort decreasing
        insertion_sort(*values, *M);

        // clean 0
        clean_zeros(values, M);
}
