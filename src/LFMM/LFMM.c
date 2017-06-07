/*
   LFMM, file: LFMM.c
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
#include <time.h>
#include <math.h>
#include "../matrix/rand.h"
#include "../matrix/matrix.h"
#include "../matrix/normalize.h"
#include "../matrix/inverse.h"
#include "../io/io_data_double.h"
#include "../io/io_data_float.h"
#include "../io/io_tools.h"
#include "print_lfmm.h"
#include "data_lfmm.h"
#include "beta.h"
#include "lfmm_k0.h"
#include "U.h"
#include "V.h"
#include "error_lfmm.h"
#include "register_lfmm.h"
#include "lfmm_algo.h"

void LFMM(LFMM_param param)
{

        // Parameters initialization
        int N2;
        param->n = 0;
        param->D = 0;

        // temporary variables
        int K = param->K;
        int n = param->n;
        int L = param->L;
        int mD = param->mD;
        int D = param->D;
        int d;

        double *perc_var;       // percentage of variances

        // random initialization
        init_random(&(param->seed));

        // count the number of lines and columns
        param->L = nb_cols_lfmm(param->input_file);
        param->n = nb_lines(param->input_file, param->L);

        param->D = nb_cols_lfmm(param->cov_file);
        N2 = nb_lines(param->cov_file, param->D);

        n = param->n;
        L = param->L;
        K = param->K;
        D = param->D;

        // check the number of lines and columns
        if (N2 != param->n) {
                Rprintf
                    ("The number of individuals of %s (%d) is different from the number"
                     " of individuals of %s (%d)\n", param->input_file,
                     param->n, param->cov_file, N2);
                error(NULL);
        }

        if (param->nd && (param->nd < 1 || param->nd > param->D))
                print_error_lfmm("specific",
                                 "(-d option). d should be between 1 and D", 0);

        // print summary of command line
        print_summary_lfmm(param);

        // allocate data memory
        // the memory is free with free_param_lfmm in the main 
        param->U = (double *)Calloc(K * n * sizeof(double), double);
        param->V = (double *)Calloc(K * L * sizeof(double), double);
        param->alpha_U = (double *)Calloc(K * sizeof(double), double);
        param->alpha_V = (double *)Calloc(K * sizeof(double), double);
        if (param->all) {
                mD = D + 1;
        } else {
                mD = 2;
        }
        param->mD = mD;
        param->beta = (double *)Calloc(mD * L * sizeof(double), double);
        param->alpha_beta = (double *)Calloc(mD * sizeof(double), double);
        perc_var = (double *)Calloc(mD + K + 1 * sizeof(double), double);

        // read of the variable file
        param->C = (double *)Calloc(n * D * sizeof(double), double);
        read_data_double(param->cov_file, n, D, param->C);
        normalize_cov(param->C, n, D);
        Rprintf("Read variable file:\n \t%s\t\tOK.\n\n", param->cov_file);

        // read of the data file
        param->dat = (float *)Calloc(n * L * sizeof(float), float);
        read_data_float(param->input_file, n, L, param->dat);

        // check that the data matrix has no constant column
        // check_constant_column(param->dat, param->n, param->L);

        // creation of the missing data matrix
        if (param->missing_data) {
                param->I = (int *)Calloc(n * L * sizeof(int), int);
                create_I(param->dat, param->I, n, L);
                inputation_freq(param->dat, param->I, n, L);
        }
        // warnings about the variables
        if (param->all) {
                Rprintf("WARNING: You launched LFMM command line with several"
                       " variables with '-a' option."
                       " The model will be\n\tlaunched with all variables at the same time.\n\n");
        } else if (!param->nd && D > 1) {
                Rprintf("WARNING: You launched LFMM command line with several"
                       " variables. The model will be\n\tlaunched sequentially"
                       " (independently) for each variable.\n\n");
        }

        Rprintf("Read genotype file:\n \t%s\t\tOK.\n", param->input_file);

        // all covariables at the same time
        if (param->all) {
                // allocate memory
                param->zscore = (double *)Calloc(L * D * sizeof(double), double);
                param->mC = (double *)Calloc(n * mD * sizeof(double), double);

                Rprintf("\n<<<<\n\t Analyse for all variables.\n\n");
                // create mC from C
                modify_C(param->C, n, D, param->mC, param->nd, param->all);

                // run LFMM
                if (K)
                        lfmm_emcmc(param);
                else
                        lfmm_k0(param);

                // write zscore
                write_zscore_double(param->output_file, L, param->zscore,
                                    mD - 1, 1, 0, K, n, param->dev, param->DIC);
                Rprintf
                    ("\tThe execution for all variables worked without error.\n>>>>\n\n");

                // only with covariable nd
        } else if (param->nd) {
                // allocate memory
                param->zscore = (double *)Calloc(L * sizeof(double), double);
                param->mC = (double *)Calloc(n * mD * sizeof(double), double);   // (N,K)
                param->nd -= 1; // modify nd to be the index of C column 

                Rprintf("\n<<<<\n\t Analyse for variable %d\n\n", param->nd + 1);
                // create mC from C
                modify_C(param->C, n, D, param->mC, param->nd, param->all);

                // run LFMM
                if (K)
                        lfmm_emcmc(param);
                else
                        lfmm_k0(param);

                // write zscore
                write_zscore_double(param->output_file, L, param->zscore,
                                    1, 0, param->nd, K, n, param->dev,
                                    param->DIC);
                Rprintf("\tThe execution for variable %d worked without error."
                       "\n>>>>\n\n", param->nd + 1);

                // each covariable sequentially
        } else {
                // allocate memory
                param->zscore = (double *)Calloc(L * sizeof(double), double);
                param->mC = (double *)Calloc(n * mD * sizeof(double), double);   // (N,K)
                // for each variable
                for (d = 0; d < param->D; d++) {
                        Rprintf("\n<<<<\n\t Analyse for variable %d\n\n", d + 1);
                        // create mC from C
                        modify_C(param->C, n, D, param->mC, d, param->all);

                        // run LFMM
                        if (K)
                                lfmm_emcmc(param);
                        else
                                lfmm_k0(param);

                        // write zscore
                        write_zscore_double(param->output_file, L,
                                            param->zscore, 1, 0, d, K, n,
                                            param->dev, param->DIC);
                        Rprintf
                            ("\tThe execution for variable %d worked without error."
                             "\n>>>>\n\n", d + 1);
                }
        }

        // free memory
        Free(perc_var);
}
