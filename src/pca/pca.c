/*
    pca, file: pca.c
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

#include "../io/io_tools.h"
#include "../io/io_data_double.h"
#include "../matrix/rand.h"
#include "calc_cov.h"
#include "print_pca.h"
#include "../matrix/normalize.h"
#include "../matrix/diagonalize.h"

// pca

void pca(char* input_file, char *output_eva_file, char *output_eve_file, 
	char *output_sdev_file, char *output_x_file,
	int *n, int *L, int *K, int c, int s)
{
        double *data;
	double *cov, *val, *vect;
	int N, M, tmp;

        // number of lines and columns
        M = nb_cols_lfmm(input_file);
        N = nb_lines(input_file, M);

	*L = M;
	*n = N;
	// correct K
	if (N < M)
		tmp = N;
	else 
		tmp = M;
	if (!(*K) || *K > tmp)
		*K = tmp;

        // print command line summary
	print_summary_pca(N, M, *K, c, s, input_file, output_eva_file, 
		output_eve_file, output_sdev_file, output_x_file);

	// allocate memory 
	data = (double *) Calloc(N * M *  sizeof(double), double);
	cov = (double *) Calloc(N * N *  sizeof(double), double);
	val = (double *) Calloc(N *  sizeof(double), double);
	vect = (double *) Calloc(N * (*K) *  sizeof(double), double);
	
	// read input_file
	read_data_double(input_file, N, M, data);

	// scale
	if (s)
		normalize_cov_I(data, N, M);
	else if (c)
		normalize_mean_I(data, N, M);

	// calculate covariance matrix
	calc_cov(data, N, M, cov); 
#ifdef USING_R
	// tout est dans le titre de la fonction,
	// check si l'utilisateur a essayé d'interrompre le programme 
	R_CheckUserInterrupt();
#endif
	// calculate eva and eve
	diagonalize(cov, N, *K, val, vect);

	// write output
	write_data_double(output_eva_file, N, 1, val);
	write_data_double(output_eve_file, N, *K, vect);

	// calculate sdev
	calc_sdev(val, N);
	write_data_double(output_sdev_file, N, 1, val);

	// calculate x into vect
	calc_x(vect, val, N, *K);
	write_data_double(output_x_file, N, *K, vect);

	// free memory
	Free(data);
	Free(cov);
	Free(val);
	Free(vect);
}
