/*
    LFMM, file: thread_var.c
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
#include "thread_var.h"
#include "LFMM.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

// thrd_var

void thrd_var(LFMM_param param, LFMM_GS_param GS_param,
              void (*fct) (), double *res, double *res2)
{
        pthread_t *thread;      // pointer to a group of threads
        int i;
        int num_thrd = param->num_thrd;

        thread = (pthread_t *) Calloc(param->num_thrd *  sizeof(pthread_t), pthread_t);
        Multithreading_lfmm_var *Ma =
            (Multithreading_lfmm_var *) malloc(param->num_thrd *
                                               sizeof(Multithreading_lfmm_var));

        /* this for loop not entered if threadd number is specified as 1 */
        for (i = 1; i < param->num_thrd; i++) {
                Ma[i] =
                    (Multithreading_lfmm_var) malloc(1 *
                                                     sizeof
                                                     (multithreading_lfmm_var));
                Ma[i]->R = param->dat;
                Ma[i]->U = param->U;
                Ma[i]->V = param->V;
                Ma[i]->C = param->mC;
                Ma[i]->beta = param->beta;
                Ma[i]->D = param->mD;
                Ma[i]->K = param->K;
                Ma[i]->N = param->n;
                Ma[i]->M = param->L;
                Ma[i]->res = 0.0;
                Ma[i]->res2 = 0.0;
                Ma[i]->num_thrd = param->num_thrd;
                Ma[i]->slice = i;
                /* creates each thread working on its own slice of i */
                if (pthread_create
                    (&thread[i], NULL, (void *)fct, (void *)Ma[i])) {
                        perror("Can't create thread");
                        Free(thread);
                        error(NULL);
                }
        }

        /* main thread works on slice 0
         *          so everybody is busy
         *                   main thread does everything if threadd number is specified as 1*/
        Ma[0] =
            (Multithreading_lfmm_var) malloc(1 *
                                             sizeof(multithreading_lfmm_var));
        Ma[0]->R = param->dat;
        Ma[0]->U = param->U;
        Ma[0]->V = param->V;
        Ma[0]->C = param->mC;
        Ma[0]->beta = param->beta;
        Ma[0]->D = param->mD;
        Ma[0]->K = param->K;
        Ma[0]->N = param->n;
        Ma[0]->M = param->L;
        Ma[0]->res = 0.0;
        Ma[0]->res2 = 0.0;
        Ma[0]->num_thrd = param->num_thrd;
        Ma[0]->slice = 0;
        /* creates each thread working on its own slice of i */
        fct(Ma[0]);

        /*main thead waiting for other thread to complete */
        for (i = 1; i < num_thrd; i++)
                pthread_join(thread[i], NULL);

        *res = 0.0;
        for (i = 0; i < num_thrd; i++)
                *res += Ma[i]->res;

        if (res2) {
                *res2 = 0.0;
                for (i = 0; i < num_thrd; i++)
                        *res2 += Ma[i]->res2;
        }

        for (i = 0; i < num_thrd; i++)
                Free(Ma[i]);
        Free(Ma);
        Free(thread);
}

// slice_mean

void slice_mean(void *G)
{
        Multithreading_lfmm_var Ma = (Multithreading_lfmm_var) G;
        double *C = Ma->C;
        double *U = Ma->U;
        double *V = Ma->V;
        double *beta = Ma->beta;
        float *R = Ma->R;
        int M = Ma->M;
        int K = Ma->K;
        int D = Ma->D;
        int N = Ma->N;
        int nb_data = N;
        int s = Ma->slice;
        int num_thrd = Ma->num_thrd;
        int from = (s * nb_data) / num_thrd;    // note that this 'slicing' works fine
        int to = ((s + 1) * nb_data) / num_thrd;        // even if SIZE is not divisible by num_thrd
        int i, j, k, d;
        double tmp1, tmp2, mean;

        mean = 0;
        for (i = from; i < to; i++) {
                for (j = 0; j < M; j++) {
                        // tmp1 = C beta
                        tmp1 = 0;
                        for (d = 0; d < D; d++)
                                tmp1 += C[i * D + d] * beta[d * M + j];
                        // tmp2 = U V
                        tmp2 = 0;
                        for (k = 0; k < K; k++)
                                tmp2 += U[k * N + i] * V[k * M + j];
                        // R - U V - C beta
                        mean += (double)(R[i * M + j]) - tmp1 - tmp2;
                }
        }
        Ma->res = mean;
}

void slice_var(void *G)
{
        Multithreading_lfmm_var Ma = (Multithreading_lfmm_var) G;
        double *C = Ma->C;
        double *U = Ma->U;
        double *V = Ma->V;
        double *beta = Ma->beta;
        float *R = Ma->R;
        int M = Ma->M;
        int K = Ma->K;
        int D = Ma->D;
        int N = Ma->N;
        int nb_data = N;
        int s = Ma->slice;
        int num_thrd = Ma->num_thrd;
        int from = (s * nb_data) / num_thrd;    // note that this 'slicing' works fine
        int to = ((s + 1) * nb_data) / num_thrd;        // even if SIZE is not divisible by num_thrd
        int i, j, k, d;
        double tmp1, tmp2, var, var2, tmp;

        var = 0.0;
        var2 = 0.0;
        for (i = from; i < to; i++) {
                for (j = 0; j < M; j++) {
                        // tmp1 = C B
                        tmp1 = 0.0;
                        for (d = 0; d < D; d++)
                                tmp1 += C[d * N + i] * beta[d * M + j];
                        // tmp2 = U V
                        tmp2 = 0.0;
                        for (k = 0; k < K; k++)
                                tmp2 += U[k * N + i] * V[k * M + j];
                        // tmp = R - UV - C B
                        tmp = ((double)(R[i * M + j]) - tmp1 - tmp2);
                        var += tmp;
                        var2 += tmp * tmp;
                }
        }
        Ma->res = var;
        Ma->res2 = var2;
}

#endif
