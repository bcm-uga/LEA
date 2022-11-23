/*
    LFMM, file: sylvester.c
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
#include <R_ext/BLAS.h>

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "sylvester.h"
#include "matrix.h"



// sylvester

void sylvester(double *A, double *B, double *C, double *X, int M, int N)
{
        double *QA = (double *)Calloc(M * M * sizeof(double), double);
        double *QB = (double *)Calloc(N * N * sizeof(double), double);
        double *Ct = (double *)Calloc(M * N * sizeof(double), double);

        int m1, m2, n1, n2;
        long int isgn = 1;
        double scale = 1.0;
        long int info;
        char trana = 'N';
        char tranb = 'N';
        long int n = (long int)N;
        long int m = (long int)M;

        schur(A, QA, M);
        schur(B, QB, N);

        // cw means colum wise/major
        // rw means row wise/major
        // Ct (cw,MxN) = QA^t(cw,MxM) * C(rw,MxN)
        for (m1 = 0; m1 < M; m1++)
                for (n1 = 0; n1 < N; n1++)
                        for (m2 = 0; m2 < M; m2++)
                                Ct[m1 + n1 * M] +=
                                    QA[m2 + m1 * M] * C[m2 * N + n1];

        // C (cw,MxN) = Ct(cw,MxN) * QB(cw,NxN)
        for (m1 = 0; m1 < M; m1++) {
                for (n1 = 0; n1 < N; n1++) {
                        C[m1 + n1 * M] = 0.0;
                        for (n2 = 0; n2 < N; n2++)
                                C[m1 + n1 * M] +=
                                    Ct[m1 + n2 * M] * QB[n2 + n1 * N];
                }
        }

        dtrsyl_(&trana, &tranb, (int *) (&isgn), (int *) (&m),
                (int *) (&n), (double *) A, (int *) (&m),
                (double *) B, (int *) (&n), (double *) C,
                (int *) (&m), (double *) (&scale), (int *) (&info) FCLEN FCLEN);

        // Ct (cw,MxN) = QA (cw, MxM) * C (cw, MxN) 
        for (m1 = 0; m1 < M; m1++) {
                for (n1 = 0; n1 < N; n1++) {
                        Ct[m1 + n1 * M] = 0.0;
                        for (m2 = 0; m2 < M; m2++)
                                Ct[m1 + n1 * M] +=
                                    QA[m1 + m2 * M] * C[m2 + n1 * M];
                }
        }

        // X (rw,MxN) = Ct (cw, MxM) * QB' (cw, MxN) 
        for (m1 = 0; m1 < M; m1++) {
                for (n1 = 0; n1 < N; n1++) {
                        X[m1 * N + n1] = 0.0;
                        for (n2 = 0; n2 < N; n2++)
                                X[m1 * N + n1] +=
                                    Ct[m1 + n2 * M] * QB[n1 + n2 * N];
                }
        }

        Free(QA);
        Free(QB);
        Free(Ct);
}

// schur

void schur(double *A, double *Q, int M)
{
        double *wr = (double *)Calloc(M * sizeof(double), double);
        double *wi = (double *)Calloc(M * sizeof(double), double);
        double *work = (double *)Calloc(3 * M * sizeof(double), double);
        long int lwork = 3 * M;
        long int info;
        long int sdim = 0;
        char jobsv = 'V';
        char sort = 'N';
        long int n = (long int)M;

        transpose_double(A, M, M);
        dgees_(&jobsv, &sort, 0, (int *) (&n), (double *) A,
               (int *) (&n), (int *) (&sdim), (double *) wr,
               (double *) wi, (double *) Q, (int *) (&n),
               (double *) work, (int *) (&lwork), 0,
               (int *) (&info) FCLEN FCLEN);

        Free(wr);
        Free(wi);
        Free(work);
}
