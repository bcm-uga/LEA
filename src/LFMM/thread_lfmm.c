/*
    LFMM, file: thread_lfmm.c
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
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

// thread_fct_lfmm

void thread_fct_lfmm(float *R, double *A, double *B, double *C, double *m,
                     double *inv_cov, double *L, int J, int K, int N, int M,
                     double *alpha, double alpha_R, int num_thrd, int mode,
                     void (*fct) ())
{
        pthread_t *thread;      // pointer to a group of threads
        int i;

        thread = (pthread_t *) Calloc(num_thrd *  sizeof(pthread_t), pthread_t);
        Multithreading_lfmm *Ma =
            (Multithreading_lfmm *) malloc(num_thrd *
                                           sizeof(Multithreading_lfmm));

        /* this for loop not entered if threadd number is specified as 1 */
        for (i = 1; i < num_thrd; i++) {
                Ma[i] =
                    (Multithreading_lfmm) malloc(1 *
                                                 sizeof(multithreading_lfmm));
                Ma[i]->R = R;
                Ma[i]->A = A;
                Ma[i]->B = B;
                Ma[i]->C = C;
                Ma[i]->m = m;
                Ma[i]->inv_cov = inv_cov;
                Ma[i]->L = L;
                Ma[i]->J = J;
                Ma[i]->K = K;
                Ma[i]->N = N;
                Ma[i]->M = M;
                Ma[i]->mode = mode;
                Ma[i]->alpha = alpha;
                Ma[i]->alpha_R = alpha_R;
                Ma[i]->num_thrd = num_thrd;
                Ma[i]->slice = i;
                /* creates each thread working on its own slice of i */
                if (pthread_create
                    (&thread[i], NULL, (void *)fct, (void *)Ma[i])) {
                        perror("Can't create thread");
                        Free(thread);
                        error(NULL);
                }
        }

        /* main thread works on slice 0 so everybody is busy
         * main thread does everything if threadd number is specified as 1*/
        Ma[0] = (Multithreading_lfmm) Calloc(1 *  sizeof(multithreading_lfmm), multithreading_lfmm);
        Ma[0]->R = R;
        Ma[0]->A = A;
        Ma[0]->B = B;
        Ma[0]->C = C;
        Ma[0]->m = m;
        Ma[0]->inv_cov = inv_cov;
        Ma[0]->L = L;
        Ma[0]->J = J;
        Ma[0]->K = K;
        Ma[0]->N = N;
        Ma[0]->M = M;
        Ma[0]->mode = mode;
        Ma[0]->alpha = alpha;
        Ma[0]->alpha_R = alpha_R;
        Ma[0]->num_thrd = num_thrd;
        Ma[0]->slice = 0;
        /* creates each thread working on its own slice of i */
        fct(Ma[0]);

        /*main thead waiting for other thread to complete */
        for (i = 1; i < num_thrd; i++)
                pthread_join(thread[i], NULL);

        for (i = 0; i < num_thrd; i++)
                Free(Ma[i]);
        Free(Ma);
        Free(thread);
}

#endif
