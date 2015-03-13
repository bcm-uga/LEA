/*
   LFMM, file: lfmm_k0.c
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

#include "print_lfmm.h"
#include "lfmm_k0.h"
#include "data_lfmm.h"
#include "print_lfmm.h"
#include "error_lfmm.h"
#include "register_lfmm.h"
#include "beta_k0.h"
#include "../matrix/matrix.h"
#include "../matrix/normalize.h"
#include "beta.h"
#include "../io/io_data_double.h"
#include <stdio.h>
#include <stdlib.h>

// lfmm_emcmc

void lfmm_k0(LFMM_param param)
{
        double *var_beta;
        double *CCt;

        // temporary parameters
        int M = param->L;
        int N = param->n;
        int D = param->mD;
        float *dat = param->dat;
        int *I = param->I;
        double *C = param->mC;
        double *zscore = param->zscore;
        double *beta = param->beta;
        int missing_data = param->missing_data;
        double perc_var;

        // allocate memory
        var_beta = (double *)Calloc(D * M *  sizeof(double), double);
        CCt = (double *)Calloc(D * D * sizeof(double), double);

        // input missing dat
        if (missing_data)
                inputation_freq(dat, I, N, M);

        // calculate C C^t
        create_CCt(CCt, C, D, N);

        // calculate E[beta] and var[beta] 
        calc_beta_k0(C, dat, beta, CCt, var_beta, M, N, D, &perc_var);

        // calculate zscore
        zscore_calc_k0(zscore, beta, var_beta, D, M);

        // check zscore
        if (check_mat(zscore, M, 0, 1))
                print_error_global("nan", NULL, 0);

        // save ED and DIC
        Rprintf("\tED: NA\t DIC: NA \n\n");

        // free memory
        Free(var_beta);
        Free(CCt);
}
