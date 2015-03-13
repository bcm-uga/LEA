/*
    calc_cov, file: calc_cov.c
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

// calc_cov

void calc_cov(double *data,  int N, int M, double *vect)
{
	int nb;
	int i1, i2, j;
	double tmp;

	for (i1 = 0; i1 < N; i1++) {
		for (i2 = 0; i2 < i1; i2++) {
			// for each couple of individuals
			tmp = 0.0;
			nb = 0;
			// calc covariance
			for (j = 0; j < M; j++) {
				if(fabs(data[i1 * M + j]) != 9 
					&& fabs(data[i2 * M + j]) != 9) {
					tmp += data[i1 * M + j] * data[i2 * M + j];
					nb++; 
				}
			}
			// check
			if (!nb) {
				Rprintf("Error: It seems that individuals %d and"
					" %d have too many missing data.\n", i1+1, i2+1);
				error(NULL);
			}
			//tmp /= nb;
			vect[i1 * N + i2] = tmp;
			vect[i2 * N + i1] = tmp;
		}
		// compute diagonal elements
		tmp = 0.0;
		nb = 0;
		for (j = 0; j < M; j++) {
			if(fabs(data[i1 * M + j]) != 9) {
				tmp += data[i1 * M + j] * data[i1 * M + j];
				nb++;
			}
		}
		if (!nb) {
			Rprintf("Error: It seems that individuals %d has"
				" too many missing data.\n", i1+1);
			error(NULL);
		}
		//tmp /= nb;
		vect[i1 * N + i1] = tmp;
	}
}

// calc_sdev

void calc_sdev(double *val, int N)
{
	int i;

	for (i = 0; i < N; i++)
		val[i] = sqrt(val[i] / (double)N);
}

// calc_x

void calc_x(double *vec, double *val, int N, int K)
{
	int i, k;

	for (i = 0; i < N; i++)
		for (k = 0; k < K; k++)
		vec[i * K + k] = vec[i * K + k] * val[k] * sqrt((double)N);
}
