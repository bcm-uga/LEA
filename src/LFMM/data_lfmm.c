/*
   LFMM, file: data_lfmm.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_lfmm.h"
#include "error_lfmm.h"
#include "register_lfmm.h"
#include "lfmm_algo.h"
#include "../matrix/cholesky.h"
#include "../matrix/inverse.h"
#include "../matrix/rand.h"
#include "../matrix/matrix.h"
#include "../io/io_tools.h"
#include "../io/io_data_double.h"
#include "../stats/gamma_distribution.h"

#ifndef WIN32
#include "thread_var.h"
#include "thread_lfmm.h"
#include "slice_lfmm.h"
#endif

// rand_matrix

void rand_matrix(double *A, double *m_A, double *inv_cov_A, double alpha_R,
                 int K, int N, int num_thrd)
{
        int i, k, kp;
        double *mu;
        double *y;
        // allocate memory
        double *L = (double *)Calloc(K * K * sizeof(double), double);

        // cholesky of inv_cov_U
        cholesky(inv_cov_A, K, L);

#ifndef WIN32
        // multithread non windows version
        if (num_thrd > 1) {
                thread_fct_lfmm(NULL, A, NULL, NULL, m_A, inv_cov_A, L,
                                0, K, N, 0, NULL, alpha_R, num_thrd, 0,
                                slice_rand);
                // uni-threaded or windows version
        } else {
#endif
                // allocate memory
                mu = (double *)Calloc(K * sizeof(double), double);
                y = (double *)Calloc(K * sizeof(double), double);
                for (i = 0; i < N; i++) {
                        for (k = 0; k < K; k++) {
                                // inv_cov_A %*% m_A
                                mu[k] = 0;
                                for (kp = 0; kp < K; kp++) {
                                        mu[k] +=
                                            inv_cov_A[k * K + kp] * m_A[kp * N +
                                                                        i];
                                }
                                // times alpha_R
                                mu[k] *= alpha_R;
                        }
                        // rand A
                        mvn_rand(mu, L, K, y);
                        for (k = 0; k < K; k++) {
                                A[k * N + i] = y[k];
                        }
                }
                // free memory
                Free(mu);
                Free(y);
#ifndef WIN32
        }
#endif
        // free memory
        Free(L);
}

// create_inv_cov

void create_inv_cov(double *inv_cov, double *alpha, double alpha_R,
                    double *A, int K, int M, int num_thrd)
{
        int k1, k2, j;
        // allocate memory
        double *tmp2 = (double *)Calloc(K * K * sizeof(double), double);

#ifndef WIN32
        // multi-threaded non windows version
        if (num_thrd > 1) {
                thread_fct_lfmm(NULL, A, NULL, NULL, NULL, tmp2, NULL,
                                0, K, 0, M, alpha, alpha_R, num_thrd, 0,
                                slice_inv_cov);
                // uni-threaded or windows version
        } else {
#endif
                // calculate cov = alphaR A %*% t(A) + diag(alpha)

                for (k1 = 0; k1 < K; k1++) {
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
                                tmp2[k1 * (K + 1)] +=
                                    (A[k1 * M + j] * A[k1 * M + j]);
                        tmp2[k1 * (K + 1)] *= alpha_R;
                        tmp2[k1 * (K + 1)] += alpha[k1];
                }
#ifndef WIN32
        }
#endif
        // inverse  
        if (K == 1) {
                inv_cov[0] = 1 / tmp2[0];
        } else {
                fast_inverse(tmp2, K, inv_cov);
        }
        // free memory
        Free(tmp2);
}

// create_m

void create_m(double *A, float *R, double *B, double *C, double *m,
              int M, int N, int J, int K, int num_thrd, int mode)
{
        int i, j, k, d;
        double *tmp_i;

        // init m with zeros
        if (!mode)
                zeros(m, K * M);

#ifndef WIN32
        // multi-threaded non windows version
        if (num_thrd >= 1) {
                thread_fct_lfmm(R, A, B, C, m, NULL, NULL,
                                J, K, N, M, NULL, 0, num_thrd, mode, slice_m);
        } else {
#endif
                // uni-threaded or windows version
                // allocate memory 
                tmp_i = (double *)Calloc(M *  sizeof(double), double);

                for (i = 0; i < N; i++) {
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
                                                m[k * N + i] +=
                                                    A[k * M + j] * tmp_i[j];
                                }
                        } else {
                                for (k = 0; k < K; k++) {
                                        for (j = 0; j < M; j++)
                                                m[k * M + j] +=
                                                    A[k * N + i] * tmp_i[j];
                                }
                        }
                }
                // free memory
                Free(tmp_i);
#ifndef WIN32
        }
#endif
}

// quantiles

void quantiles(double *dist, double *prob, int n, int p, double *res)
{
        int *index = (int *)Calloc(n * sizeof(int), int);
        int j, jm, jp;

        sort_index(dist, index, n);
        for (j = 0; j < p; j++) {
                jm = floor(n * prob[j]);
                jp = ceil(n * prob[j]);
                res[j] = (dist[index[jm]] + dist[index[jp]]) / 2;
        }

        Free(index);
}

// lambda

double lambda(double *p, int n)
{
        //  qchisq(.5, df=1)/ median(qchisq(p.values, df=1))
        // double dot5 = 0.4549364;
        double *qchisq = (double *)Calloc(41 * sizeof(double), double);
        int i;
        double *pp = (double *)Calloc(41 * sizeof(double), double);
        double *q = (double *)Calloc(41 * sizeof(double), double);
        double *dist = (double *)Calloc(n * sizeof(double), double);
        double res;

        pp[0] = 0.5;
        for (i = 1; i < 41; i++)
                pp[i] = pp[i - 1] + 0.01;

        for (i = 0; i < n; i++)
                dist[i] = p[i] * p[i];

        quantiles(dist, pp, n, 41, q);
        print_data_double(q, 1, 41);

        for (i = 0; i < 41; i++) {
                qchisq[i] = q[i] / quantile_Gamma_Distribution(pp[i], .5);
        }

        res = median(qchisq, 41);

        Free(qchisq);
        Free(pp);
        Free(q);
        Free(dist);

        return res;
}

// pvalue_qvalue

void pvalue_qvalue(double *pvalues, double *qvalues, int n)
{
        int *index = (int *)Calloc(n * sizeof(int), int);
        int i;

        // sort pvalue table
        sort_index(pvalues, index, n);

        for (i = 0; i < n; i++) {
                qvalues[index[i]] = pvalues[index[i]]
                    * (double)(n) / (double)(i + 1);
                if (qvalues[index[i]] > 1.0)
                        qvalues[index[i]] = 1.0;
        }
        Free(index);
}

// zscore_calc

void zscore_calc(double *zscore, double *sum, double *sum2, int n, int cur,
                 int D)
{
        int i;
        double var;
        double *r = (double *)Calloc(n * (D - 1) * sizeof(double), double);
        double *m = (double *)Calloc(n * (D - 1) * sizeof(double), double);

        for (i = n; i < D * n; i++) {
                // calculate var beta
                var =
                    ((sum2[i] - sum[i] * sum[i] / (double)cur) / (double)(cur -
                                                                          1));

                // calculate zscore beta
                zscore[(i - n)] = (sum[i]) / (sqrt(var) * (double)cur);
                r[i - n] = var;
                m[i - n] = sum[i] / cur;
        }

        Free(r);
        Free(m);
}

// update_sum

void update_sum(double *beta, double *sum, int n)
{
        int i;

        for (i = 0; i < n; i++)
                sum[i] += beta[i];
}

// update_sum2

void update_sum2(double *beta, double *sum2, int n)
{
        int i;

        for (i = 0; i < n; i++)
                sum2[i] += (beta[i] * beta[i]);
}

// update_deviance

void update_deviance(double *deviance, int cur, double var, double thrd_m)
{
        double mean = thrd_m / var;

        *deviance = ((double)(cur - 1) * (*deviance) + mean) / (double)cur;
}

// update_mean

void update_mean(double *beta, double *mean, int n, int cur)
{
        int i;

        for (i = 0; i < n; i++)
                mean[i] = ((double)(cur - 1) * mean[i] + beta[i]) / (double)cur;
}

// modify_C

void modify_C(double *C, int N, int nD, double *Cpp, int d, int all)
{
        int i, nd;

        if (all) {
                for (i = 0; i < N; i++)
                        Cpp[i] = 1.0;
                for (i = 0; i < N; i++)
                        for (nd = 0; nd < nD; nd++)
                                Cpp[(nd + 1) * (N) + i] = C[i * nD + nd];

        } else {
                for (i = 0; i < N; i++)
                        Cpp[i] = 1.0;
                for (i = 0; i < N; i++)
                        Cpp[N + i] = C[i * nD + d];
        }
}

// write_DIC

void write_DIC(char *file_data, double deviance, double DIC)
{
        FILE *file = NULL;

        file = fopen_write(file_data);

        fprintf(file, "Deviance DIC\n");
        fprintf(file, "%G %G\n", deviance, DIC);
        fclose(file);
}

// write_zscore_double

void write_zscore_double(char *output_file, int M, double *zscore, int D,
                         int all, int nd, int K, int N, double dev, double DIC)
{
        FILE *file = NULL;
        FILE *file_dic = NULL;
        int j;
        int d;
        char zscore_file[512];
        char dic_file[512];
        double *pvalues = (double *)Calloc(M * sizeof(double), double);
        // double* qvalues = (double *)Calloc(M * sizeof(double), double);

        if (all) {
                // DIC file
                snprintf(dic_file, 512, "%s_a.%d.dic", output_file, K);
                file_dic = fopen_write(dic_file);

                fprintf(file_dic, "K\t\t\t%d\n", K);
                fprintf(file_dic, "D\t\t\t%d\n", D);
                if (K) {
                        fprintf(file_dic, "Deviance\t\t%10.10G\n", dev);
                        fprintf(file_dic, "DIC\t\t\t%10.10G\n", DIC);
                } else
                        fprintf(file_dic, "Deviance\t\tNa\nDIC\t\t\tNa\n");

                Rprintf
                    ("\tThe statistics for the run are registered in:\n \t\t%s.\n"
                     "\n", dic_file);

                Rprintf("\t-------------------------\n");
                for (d = 0; d < D; d++) {
                        // calculate pvalues
                        for (j = 0; j < M; j++)
                                pvalues[j] =
                                    (double)
                                    zscore2pvalue_student(fabs
                                                          (zscore[d * M + j]),
                                                          N - D);
                        // calculate qvalues
                        // pvalue_qvalue(pvalues, qvalues, M);

                        // calculate lambda
                        /*
                           l = lambda(pvalues, M);
                           Rprintf("\tLambda for variable %d:\t%3.5G\n\n", d+1, l);
                           // write in file
                           fprintf(file_dic, "lambda_v%d\t\t%3.5G\n", d + 1, l);
                         */
                        // write file name 
                        snprintf(zscore_file, 512, "%s_a%d.%d.zscore",
                                 output_file, d + 1, K);
                        // and write 
                        file = fopen_write(zscore_file);
                        for (j = 0; j < M; j++) {
                                fprintf(file, "%3.6G %3.6G %3.6G",
                                        (double)zscore[d * M + j],
                                        (double)(-log10(pvalues[j])),
                                        pvalues[j]);
                                //                              (double)(-log10(qvalues[j])), qvalues[j]);
                                fprintf(file, "\n");
                        }
                        fclose(file);
                        Rprintf
                            ("\tThe zscores for variable %d are registered in:\n \t\t%s.\n"
                             "\tThe columns are: zscores, -log10(p-values), p-values.\n"
                             "\n", d + 1, zscore_file);
                        Rprintf("\t-------------------------\n");
                }
                // write percentage
                /*
                   for (d = 0; d < D+1; d++)
                   fprintf(file_dic, "percentage_var%d\t\t%3.3G\n", d, prec_var[d+1]);
                   for (d = 0; d < K; d++)
                   fprintf(file_dic, "percentage_factor%d\t\t%3.3G\n", d+1, prec_var[d+2+D]);
                   fprintf(file_dic, "percentage_residual\t%3.3G\n", prec_var[0]);
                 */
                fclose(file_dic);
        } else {
                // DIC file
                snprintf(dic_file, 512, "%s_s%d.%d.dic", output_file, nd + 1,
                         K);
                file_dic = fopen_write(dic_file);
                fprintf(file_dic, "K\t\t\t%d\n", K);
                fprintf(file_dic, "D\t\t\t%d\n", D);
                if (K) {
                        fprintf(file_dic, "Deviance\t\t%10.10G\n", dev);
                        fprintf(file_dic, "DIC\t\t\t%10.10G\n", DIC);
                } else
                        fprintf(file_dic, "Deviance\t\tNa\nDIC\t\t\tNa\n");

                // calculate pvalues
                for (j = 0; j < M; j++)
                        pvalues[j] =
                            (double)zscore2pvalue_student(fabs(zscore[j]),
                                                          N - D);
                // calculate qvalues
                //      pvalue_qvalue(pvalues, qvalues, M);

                Rprintf
                    ("\tThe statistics for the run are registered in:\n \t\t%s.\n"
                     "\n", dic_file);

                // calculate lambda
                /*
                   l = lambda(pvalues, M);
                   Rprintf("\t-------------------------\n");
                   Rprintf("\tLambda for variable %d:\t%3.5G\n\n", nd + 1, l);
                   fprintf(file_dic, "lambda_v%d\t\t%3.5G\n", nd + 1, l);
                 */
                // write percentage
                /*
                   for (d = 0; d < D+1; d++)
                   fprintf(file_dic, "percentage_var%d\t\t%3.3G\n", d, prec_var[d+1]);
                   for (d = 0; d < K; d++)
                   fprintf(file_dic, "percentage_factor%d\t\t%3.3G\n", d+1, prec_var[d+2+D]);
                   fprintf(file_dic, "percentage_residual\t%3.3G\n", prec_var[0]);
                 */
                fclose(file_dic);

                // write file name
                snprintf(zscore_file, 512, "%s_s%d.%d.zscore", output_file,
                         nd + 1, K);
                // and write it
                file = fopen_write(zscore_file);
                for (j = 0; j < M; j++) {
                        fprintf(file, "%3.6G %3.6G %3.6G", (double)zscore[j],
                                (double)(-log10(pvalues[j])), pvalues[j]);
                        //                      (double)(-log10(qvalues[j])), qvalues[j]);
                        fprintf(file, "\n");
                }
                fclose(file);
                Rprintf
                    ("\tThe zscores for variable %d are registered in:\n \t\t%s.\n"
                     "\tThe columns are: zscores, -log10(p-values), p-values.\n"
                     "\n", nd + 1, zscore_file);
                Rprintf("\t-------------------------\n");
        }

        Free(pvalues);
        // Free(qvalues);
}

// var_data

double var_data(LFMM_param param, LFMM_GS_param GS_param)
{
        double mean, mean2, tmp1, tmp2, tmp;
        int i, j, d, k;
        int N = param->n;
        int M = param->L;
        int K = param->K;
        int D = param->mD;
        /*
           thrd_var(R,U,V,C,beta,K,D,M,N,num_thrd,slice_mean,0,&mean,0);
           mean /= N*M;
         */
#ifndef WIN32
        if (param->num_thrd > 1) {
                thrd_var(param, GS_param, slice_var, &mean, &mean2);
        } else {
#endif
                mean = 0.0;
                mean2 = 0.0;
                for (i = 0; i < N; i++) {
                        for (j = 0; j < M; j++) {
                                tmp1 = 0.0;
                                for (d = 0; d < D; d++)
                                        tmp1 +=
                                            param->mC[d * N +
                                                      i] * param->beta[d * M +
                                                                       j];
                                tmp2 = 0.0;
                                for (k = 0; k < K; k++)
                                        tmp2 +=
                                            param->U[k * N +
                                                     i] * param->V[k * M + j];
                                tmp =
                                    ((double)(param->dat[i * M + j]) - tmp1 -
                                     tmp2);
                                mean += tmp;    //(double)(R[i*M+j])-tmp1-tmp2)
                                mean2 += tmp * tmp;     //((double)(R[i*M+j])-tmp1-tmp2)*
                                // ((double)(R[i*M+j])-tmp1-tmp2);
                        }
                }
#ifndef WIN32
        }
#endif

        GS_param->thrd_m2 = mean2;

        return (mean2 - mean * mean / (N * M)) / (N * M - 1);
}

// var_data_inputation
/*
double var_data_inputation(float *R, int *I, double *U, double *V, double *C, 
	double *beta, int N, int M, int K, int D, double *thrd_m2, int num_thrd)
{
	double mean, mean2, tmp1, tmp2, tmp;
	int i, j, d, k;

# ifndef WIN32
	if (num_thrd > 1) {
		thrd_var(R, U, V, C, beta, K, D, M, N, num_thrd, slice_var, 0, &mean,
				&mean2);
	} else {
# endif
		mean = 0.0;
		mean2 = 0.0;
		for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
				tmp1 = 0.0;
				for (d = 0; d < D; d++)
					tmp1 += C[i * D + d] * beta[d * M + j];
				tmp2 = 0.0;
				for (k = 0; k < K; k++)
					tmp2 += U[k * N + i] * V[k * M + j];
				tmp = ((double)(R[i * M + j]) - tmp1 - tmp2);
				if (!I[i * M + j])	
					R[i * M + j] = (float)(tmp1 + tmp2);
				mean += tmp;     //(double)(R[i*M+j])-tmp1-tmp2)
				mean2 += tmp * tmp;  	//((double)(R[i*M+j])-tmp1-tmp2)*
				// ((double)(R[i*M+j])-tmp1-tmp2);
			}
		}
# ifndef WIN32
	}
# endif

	*thrd_m2 = mean2;

	return (mean2 - mean * mean / (N * M)) / (N * M - 1);
}
*/

// inputation

void inputation_lfmm(float *R, double *U, double *V, double *C, double *beta,
                     int *I, int N, int M, int K, int D)
{
        int i, j, k, d;
        double tmp1, tmp2;

        // for each data in R
        for (i = 0; i < N; i++) {
                for (j = 0; j < M; j++) {
                        // if missing data
                        if (!I[i * M + j]) {
                                // U V + C beta
                                tmp1 = 0.0;
                                for (d = 0; d < D; d++)
                                        tmp1 += C[i * D + d] * beta[d * M + j];
                                tmp2 = 0.0;
                                for (k = 0; k < K; k++)
                                        tmp2 += U[k * N + i] * V[k * M + j];
                                // imput
                                R[i * M + j] = (float)(tmp1 + tmp2);
                        }
                }
        }
}

// inputation_freq

void inputation_freq(float *R, int *I, int N, int M)
{
        int i, j, nb;
        double freq;

        // for each data in R
        for (j = 0; j < M; j++) {
                freq = 0.0;
                nb = 0;
                for (i = 0; i < N; i++) {
                        if (I[i * M + j]) {
                                freq += R[i * M + j];
                                nb++;
                        }
                }
                if (nb) 
                    freq /= nb;
                for (i = 0; i < N; i++) {
                        if (!I[i * M + j]) {
                                R[i * M + j] = rand_binary(freq)
                                    + rand_binary(freq);
                        }
                }
        }
}
