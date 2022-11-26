/*
    matrix, file: inverse.c
    Copyright (C) 2012 Eric Frichot / 2022 Olivier Francois

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
#include "inverse.h"
#include "matrix.h"





void fast_inverse(double *A, int D, double *inv)
{
        long int *pivot = (long int *)Calloc(D + 1 * sizeof(long int), long int);
        double *tmp = (double *)Calloc(D * D * sizeof(double), double);
        long int info;
        long int Dp = (long int)D;
        long int size = Dp * Dp;

        copy_vect(A, inv, D * D);

        dgetrf_((int *) (&Dp), (int *) (&Dp), (double *) inv,
                (int *) (&Dp), (int *) pivot, 
                (int *) (&info) FCONE);
        dgetri_((int *) (&Dp), (double *) inv, (int *) (&Dp),
                (int *) pivot, (double *) tmp, (int *) (&size),
                (int *) (&info));

        Free(tmp);
        Free(pivot);
}

// detrm (from the web)

double detrm(double *a, int k)
{
        double s = 1, det = 0;
        int i, j, m, n, c;
        double *b = (double *)Calloc((k - 1) * (k - 1) * sizeof(double), double);

        if (k == 1) {
                Free(b);
                return (a[0]);
        } else {
                det = 0;
                for (c = 0; c < k; c++) {
                        m = 0;
                        n = 0;
                        for (i = 0; i < k; i++) {
                                for (j = 0; j < k; j++) {
                                        // b[ i*k + j ] = 0;
                                        if (i != 0 && j != c) {
                                                b[m * (k - 1) + n] =
                                                    a[i * k + j];
                                                if (n < (k - 2))
                                                        n++;
                                                else {
                                                        n = 0;
                                                        m++;
                                                }
                                        }
                                }
                        }
                        det = det + s * (a[c] * detrm(b, k - 1));
                        s = -1 * s;
                }
        }

        Free(b);
        return (det);
}

// cofact (from the web)

void cofact(double *num, int f, double *inv)
{
        double *b = (double *)Calloc((f - 1) * (f - 1) * sizeof(double), double);
        double *fac = (double *)Calloc(f * f * sizeof(double), double);
        int p, q, m, n, i, j;

        for (q = 0; q < f; q++) {
                for (p = 0; p < f; p++) {
                        m = 0;
                        n = 0;
                        for (i = 0; i < f; i++) {
                                for (j = 0; j < f; j++) {
                                        // b[ i*(f-1)+j ] = 0;
                                        if (i != q && j != p) {
                                                b[m * (f - 1) + n] =
                                                    num[i * f + j];
                                                if (n < (f - 2))
                                                        n++;
                                                else {
                                                        n = 0;
                                                        m++;
                                                }
                                        }
                                }
                        }

                        fac[q * f + p] = pow(-1, q + p) * detrm(b, f - 1);
                }
        }

        trans(num, fac, f, inv);
        Free(b);
        Free(fac);
}

// trans (from the web)

void trans(double *num, double *fac, int r, double *inv)
{
        int i, j;
        double d;
        double *b = (double *)Calloc(r * r * sizeof(double), double);

        for (i = 0; i < r; i++) {
                for (j = 0; j < r; j++) {
                        b[i * r + j] = fac[j * r + i];
                }
        }

        d = detrm(num, r);

        for (i = 0; i < r; i++) {
                for (j = 0; j < r; j++) {
                        inv[i * r + j] = b[i * r + j] / d;
                }
        }

        Free(b);
}
