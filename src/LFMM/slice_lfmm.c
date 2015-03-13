/*
   LFMM, file: slice_lfmm.c
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

#ifndef WIN32

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "thread_lfmm.h"
#include "../matrix/rand.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

// slice_m

void slice_m(void *G)
{
        Multithreading_lfmm Ma = (Multithreading_lfmm) G;
        float *R = Ma->R;
        double *A = Ma->A;
        double *B = Ma->B;
        double *C = Ma->C;
        double *m = Ma->m;
        int N = Ma->N;
        int M = Ma->M;
        int K = Ma->K;
        int J = Ma->J;
        int mode = Ma->mode;
        int nb_data, s, num_thrd, from, to;
        int i, j, k, d;
        double *tmp_i;
        double *tmp_j;
        /*
           int nb_data = N;
           int s = Ma->slice;
           int num_thrd = Ma->num_thrd;
           int from = (s * nb_data) / num_thrd; // note that this 'slicing' works fine
           int to = ((s + 1) * nb_data) / num_thrd;     // even if SIZE is not divisible by num_thrd
           int i, j, k, d;

           // allocate memory
           double *tmp_i = (double *) Calloc(M *  sizeof(double), double);
         */

        /*
           for (i = from; i < to; i++) {
           // calculate tmp_i = R - B'C
           for (j = 0; j < M; j++)
           tmp_i[j] = (double)(R[i * M + j]);
           for (d = 0; d < J; d++) {
           for (j = 0; j < M; j++)
           tmp_i[j] -= B[d * N + i] * C[d * M + j];
           }
           // calculate tmp_i * A'
           if (mode) {
           for (k = 0; k < K; k++) {
           m[k * N + i] = 0;
           for (j = 0; j < M; j++)
           m[k * N + i] += A[k * M + j] * tmp_i[j];
           }
           } else {
           for (k = 0; k < K; k++) {
           for (j = 0; j < M; j++)
           m[k * M + j] += A[k * N + i] * tmp_i[j];
           }
           }
           } */
        if (mode) {
                nb_data = N;
                s = Ma->slice;
                num_thrd = Ma->num_thrd;
                from = (s * nb_data) / num_thrd;        // note that this 'slicing' works fine
                to = ((s + 1) * nb_data) / num_thrd;    // even if SIZE is not divisible by num_thrd

                tmp_i = (double *)Calloc(M *  sizeof(double), double);

                for (i = from; i < to; i++) {
                        // calculate tmp_i = R - B'C
                        for (j = 0; j < M; j++)
                                tmp_i[j] = (double)(R[i * M + j]);
                        for (d = 0; d < J; d++) {
                                for (j = 0; j < M; j++)
                                        tmp_i[j] -= B[d * N + i] * C[d * M + j];
                        }
                        // calculate tmp_i * A'
                        for (k = 0; k < K; k++) {
                                m[k * N + i] = 0;
                                for (j = 0; j < M; j++)
                                        m[k * N + i] += A[k * M + j] * tmp_i[j];
                        }
                }

                // free memory
                Free(tmp_i);

        } else {
                nb_data = M;
                s = Ma->slice;
                num_thrd = Ma->num_thrd;
                from = (s * nb_data) / num_thrd;        // note that this 'slicing' works fine
                to = ((s + 1) * nb_data) / num_thrd;    // even if SIZE is not divisible by num_thrd

                tmp_j = (double *)Calloc(N *  sizeof(double), double);

                for (j = from; j < to; j++) {
                        // calculate tmp_i = R - B'C
                        for (i = 0; i < N; i++)
                                tmp_j[i] = (double)(R[i * M + j]);
                        for (d = 0; d < J; d++) {
                                for (i = 0; i < N; i++)
                                        tmp_j[i] -= B[d * N + i] * C[d * M + j];
                        }
                        // calculate tmp_i * A'
                        for (k = 0; k < K; k++) {
                                m[k * M + j] = 0;
                                for (i = 0; i < N; i++)
                                        m[k * M + j] += A[k * N + i] * tmp_j[i];
                        }
                }

                // free memory
                Free(tmp_j);
        }

}

// slice_rand

void slice_rand(void *G)
{
        Multithreading_lfmm Ma = (Multithreading_lfmm) G;
        double *m = Ma->m;
        double *inv_cov = Ma->inv_cov;
        double *A = Ma->A;
        double *L = Ma->L;
        double alpha_R = Ma->alpha_R;
        int N = Ma->N;
        int K = Ma->K;
        int nb_data = N;
        int s = Ma->slice;
        int num_thrd = Ma->num_thrd;
        int from = (s * nb_data) / num_thrd;    // note that this 'slicing' works fine
        int to = ((s + 1) * nb_data) / num_thrd;        // even if SIZE is not divisible by num_thrd
        int i, k, kp;

        // allocate memory
        double *mu = (double *)Calloc(K * sizeof(double), double);
        double *y = (double *)Calloc(K * sizeof(double), double);

        for (i = from; i < to; i++) {
                for (k = 0; k < K; k++) {
                        // inv_cov %*% m
                        mu[k] = 0;
                        for (kp = 0; kp < K; kp++) {
                                mu[k] += inv_cov[k * K + kp] * m[kp * N + i];
                        }
                        // times alpha_R
                        mu[k] *= alpha_R;
                }
                // rand A
                mvn_rand(mu, L, K, y);
                for (k = 0; k < K; k++)
                        A[k * N + i] = y[k];
        }
        // free memory
        Free(mu);
        Free(y);
}

// slice_inv_cov

void slice_inv_cov(void *G)
{
        Multithreading_lfmm Ma = (Multithreading_lfmm) G;
        double *tmp2 = Ma->inv_cov;
        double *A = Ma->A;
        int M = Ma->M;
        int K = Ma->K;
        double *alpha = Ma->alpha;
        double alpha_R = Ma->alpha_R;
        int nb_data = K;
        int s = Ma->slice;
        int num_thrd = Ma->num_thrd;
        int from = (s * nb_data) / num_thrd;    // note that this 'slicing' works fine
        int to = ((s + 1) * nb_data) / num_thrd;        // even if SIZE is not divisible by num_thrd
        int j, k1, k2;

        // calculate cov_U (tmp2) = alphaR A %*% t(A) + diag(alpha)
        for (k1 = from; k1 < to; k1++) {
                for (k2 = 0; k2 < k1; k2++) {
                        tmp2[k1 * K + k2] = 0;
                        for (j = 0; j < M; j++)
                                tmp2[k1 * K + k2] +=
                                    (A[k1 * M + j] * A[k2 * M + j]);
                        tmp2[k1 * K + k2] *= alpha_R;
                        tmp2[k2 * K + k1] = tmp2[k1 * K + k2];
                }
                tmp2[k1 * (K + 1)] = 0;
                for (j = 0; j < M; j++)
                        tmp2[k1 * (K + 1)] += (A[k1 * M + j] * A[k1 * M + j]);
                tmp2[k1 * (K + 1)] *= alpha_R;
                tmp2[k1 * (K + 1)] += alpha[k1];
        }
}

#endif
