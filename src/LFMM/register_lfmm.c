/*
    LFMM, file: register_lfmm.c
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
#include "register_lfmm.h"
#include "print_lfmm.h"
#include "error_lfmm.h"
#include "../io/io_tools.h"

// init_param_lfmm

void init_param_lfmm(LFMM_param param)
{
        // default values 
        param->nd = 0;
        param->D = 0;
        param->mD = 0;
        param->all = 0;
        param->K = 0;
        param->burn = 5000;
        param->Niter = 10000;
        param->num_thrd = 1;
        param->missing_data = 0;
        strcpy(param->output_file, "");
        param->seed = -1;
        param->noise_epsilon = 1e-3;
        param->b_epsilon = 1000;
        param->init = 1;

}

// analyse_param_lfmm

void analyse_param_lfmm(int argc, char *argv[], LFMM_param param)
{
        // temporary variables 
        int i, g_data = 0, g_cov = 0;
        char *tmp_file;
        int g_d = 0;

        // read each command line argument
        for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
                                // all 
                        case 'a':
                                param->all = 1;
                                break;
                                // nd-th variable
                        case 'd':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "d (numerous of the covariable)",
                                                         0);
                                param->nd = atoi(argv[i]);
                                g_d = 1;
                                break;
                                // number of latent factors
                        case 'K':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "K (number of latent factors)",
                                                         0);
                                param->K = atoi(argv[i]);
                                break;
                                // number of individuals (deprecated)
                        case 'n':
                                i++;
                                if (argv[i][0] == '-')
                                        i--;
                                Rprintf
                                    ("Warning: '-n' option is not necessary, (from LFMM v1.3)."
                                     " The number of individuals is automatically computed.\n\n");
                                break;
                                // number of loci (deprecated)
                        case 'L':
                                i++;
                                if (argv[i][0] == '-')
                                        i--;
                                Rprintf
                                    ("Warning: '-L' option is not necessary, (from LFMM v1.3)."
                                     " The number of loci is automatically computed.\n\n");
                                break;
                                // number of variables (deprecated)
                        case 'D':
                                i++;
                                if (argv[i][0] == '-')
                                        i--;
                                Rprintf
                                    ("Warning: '-D' option is not necessary (from LFMM v1.3)."
                                     " The number of environmental variables is automatically computed.\n\n");
                                break;
                                // number of iterations in the GS
                        case 'i':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "i (number of iterations in the GS)",
                                                         0);
                                param->Niter = atoi(argv[i]);
                                break;
                                // seed
                        case 's':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "s (seed number)", 0);
                                param->seed = atoll(argv[i]);
                                break;
                                // missing data
                        case 'm':      // global
                                param->missing_data = 1;
                                break;
                                // help
                        case 'h':      // global
                                print_help_lfmm();
                                error(NULL);
                                break;
                                // licence
                        case 'l':      // global
                                print_licence_lfmm();
                                error(NULL);
                                break;
                                // burnin
                        case 'b':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "b (burn parameter in the GS)",
                                                         0);
                                param->burn = atoi(argv[i]);
                                break;
                                // CPU
                        case 'p':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "p (number of processes used)",
                                                         0);
                                param->num_thrd = atoi(argv[i]);
                                break;
                                // output 
                        case 'o':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "o (output file with z-scores)",
                                                         0);
                                strcpy(param->output_file, argv[i]);
                                break;
                                // genotypes 
                        case 'x':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "x (genotype file)",
                                                         0);
                                g_data = 1;
                                strcpy(param->input_file, argv[i]);
                                break;
                                // genotypes (previous version (maintained)
                        case 'g':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "x (genotype file)",
                                                         0);
                                g_data = 1;
                                strcpy(param->input_file, argv[i]);
                                break;
                                // variables
                        case 'v':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_lfmm("cmd",
                                                         "v (variable file)",
                                                         0);
                                g_cov = 1;
                                strcpy(param->cov_file, argv[i]);
                                break;
                        default:
                                print_error_lfmm("basic", NULL, 0);
                        }
                } else {
                        print_error_lfmm("basic", NULL, 0);
                }
        }

        // checks 
        if (!g_data)
                print_error_lfmm("option", "-x genotype_file", 0);

        if (!g_cov)
                print_error_lfmm("option", "-v variable_file", 0);

        if (param->all && param->nd) {
                print_error_lfmm("specific",
                                 "-a (to run LFMM on all covariables at the same time) and -d (to run LFMM"
                                 " on a specific variable) cannot be provided in the same command line.",
                                 0);
        }

        if (g_d && param->nd <= 0)
                print_error_lfmm("missing", NULL, 0);

        if (param->K < 0 || param->num_thrd <= 0 || param->burn <= 0
            || param->Niter <= 0)
                print_error_lfmm("missing", NULL, 0);

        if (param->burn >= param->Niter) {
                print_error_lfmm("specific",
                                 "the number of iterations for burnin "
                                 "(b) is greater than the number total of iterations (i)",
                                 0);
        }
        // write output file name
        tmp_file = remove_ext(param->input_file, '.', '/');
        if (!strcmp(param->output_file, ""))
                strcpy(param->output_file, tmp_file);
        Free(tmp_file);
}

// free_param 

void free_param_lfmm(LFMM_param param)
{
        // alpha_beta
        if (param->alpha_beta)
                Free(param->alpha_beta);
        // alpha_U
        if (param->alpha_U)
                Free(param->alpha_U);
        // alpha_V
        if (param->alpha_V)
                Free(param->alpha_V);
        // I
        if (param->I)
                Free(param->I);
        // U
        if (param->U)
                Free(param->U);
        // V
        if (param->V)
                Free(param->V);
        // dat
        if (param->dat)
                Free(param->dat);
        // beta
        if (param->beta)
                Free(param->beta);
        // C
        if (param->C)
                Free(param->C);
        // mC
        if (param->mC)
                Free(param->mC);
        // zscore
        if (param->zscore)
                Free(param->zscore);
}
