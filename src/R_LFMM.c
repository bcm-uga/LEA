/*
    LFMM, file: main.c
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
#include "LFMM/LFMM.h"
#include "LFMM/register_lfmm.h"

#include "R_LFMM.h"

void R_LFMM(char **R_input_file, char **R_output_file, char **R_cov_file,
            int *R_n, int *R_L, int *R_D, int *R_nd,
            int *R_K, int *R_Niter, int *R_burn,
            int *R_num_CPU, long long *R_seed, int *R_missing_data, int *R_all,
            double *R_dic, double *R_dev, double *noise_epsilon,
            double *b_epsilon, int *init)
{
        // for random in R
        GetRNGstate();

        // parameters allocation
        lfmm_param *param = (lfmm_param *) calloc(1, sizeof(lfmm_param));

        // Parameters initialization
        init_param_lfmm(param);
        param->K = *R_K;
        param->nd = *R_nd;
        param->Niter = *R_Niter;
        param->burn = *R_burn;
        param->num_thrd = *R_num_CPU;
        param->init = *init;
        param->noise_epsilon = *noise_epsilon;
        param->b_epsilon = *b_epsilon;
        param->missing_data = *R_missing_data;
        param->seed = *R_seed;
        param->all = *R_all;
        strcpy(param->output_file, *R_output_file);
        strcpy(param->input_file, *R_input_file);
        strcpy(param->cov_file, *R_cov_file);
        param->n = *R_n;
        param->L = *R_L;

        // run
        LFMM(param);

        // output
        *R_D = param->D;
        *R_n = param->n;
        *R_L = param->L;
        *R_dic = param->DIC;
        *R_dev = param->dev;

        // free memory
        free_param_lfmm(param);
        free(param);

        // for random in R
        PutRNGstate();
}
