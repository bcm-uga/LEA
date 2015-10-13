/*
    NMF, file: main.c
    Copyright (C) 2013 François Mathieu, Eric Frichot

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "sNMF/sNMF.h"
#include "sNMF/register_snmf.h"

#include "R_sNMF.h"

void R_sNMF(char **R_genotype_file, int *R_K, double *R_alpha, double *R_tol,
            double *R_percentage, int *R_iteration, long long *R_seed,
            int *R_ploidy, int *R_num_proc, char **R_input_file_Q,
            char **R_output_file_Q, char **R_output_file_G, int *I,
            double *all_ce, double *masked_ce, int *n, int *L)
{
        // for random in R
        GetRNGstate();

        // parameters allocation
        sNMF_param param = (sNMF_param) calloc(1, sizeof(snmf_param));

        // init parameters
        init_param_snmf(param);
        param->K = *R_K;
        param->seed = *R_seed;
        param->maxiter = *R_iteration;
        param->num_thrd = *R_num_proc;
        param->tolerance = *R_tol;
        param->pourcentage = *R_percentage;
        param->alpha = *R_alpha;
        param->I = *I;
        param->m = *R_ploidy;
        // param->missing_data = R_;
        param->seed = *R_seed;
        strcpy(param->input_file, *R_genotype_file);
        strcpy(param->input_file_Q, *R_input_file_Q);
        strcpy(param->output_file_Q, *R_output_file_Q);
        strcpy(param->output_file_F, *R_output_file_G);

        // run
        sNMF(param);

        *all_ce = param->all_ce;
        *masked_ce = param->masked_ce;
        *R_seed = param->seed;
        *n = param->n;
        *L = param->L;

        // free memory
        free_param_snmf(param);
        free(param);

        // for random in R
        PutRNGstate();
}
