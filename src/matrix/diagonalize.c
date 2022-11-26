/*
    matrix, file: diagonalize.c
    Copyright (C) 2012 Eric Frichot

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
#include <R_ext/Lapack.h>

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "diagonalize.h"




#include "../io/io_data_double.h"

void diagonalize(double *cov, int N, int K, double *val, double *vect)
{
        long int n = (long int)N;
        long int M = (long int)K;
        double abstol = 1e-10;
        long int *supp = (long int *)Calloc(2 * N *  sizeof(long int), long int);
        long int lwork = 26 * N;
        double *work = (double *)Calloc(lwork * sizeof(double), double);
        long int liwork = 10 * N;
        long int *iwork = (long int *)Calloc(liwork * sizeof(double), double);
        long int info;
        double vl = 0.0, vu = 0.0;
        char jobz = 'V', range = 'I', uplo = 'U';
        long int il = (long int)N - K + 1;
        long int ul = (long int)N;
        double *valp = (double *)Calloc(N * sizeof(double), double);
        double *vectp = (double *)Calloc(N * N * sizeof(double), double);
        int i, k;

        dsyevr_((char *)(&jobz), (char *)(&range), (char *)(&uplo),
                (int *) (&n), (double *) cov, (int *) (&n),
                (double *) (&vl), (double *) (&vu), (int *) (&il),
                (int *) (&ul), (double *) (&abstol), (int *) (&M),
                (double *) valp, (double *) vectp, (int *) (&n),
                (int *) supp, (double *) work, (int *) (&lwork),
                (int *) iwork, (int *) (&liwork), 
                (int *) (&info) FCONE FCONE FCONE);

        // copy results
        for (k = 0; k < K; k++) {
                val[k] = valp[K - (k + 1)];
                if (val[k] < 0 && fabs(val[k]) < 1e-10)
                        val[k] = 0;
        }

        for (k = 0; k < K; k++)
                for (i = 0; i < N; i++)
                        vect[i * K + k] = vectp[(K - (k + 1)) * N + i];

        Free(valp);
        Free(vectp);
        Free(supp);
        Free(work);
        Free(iwork);
}
