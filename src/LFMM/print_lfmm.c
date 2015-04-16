/*
    LFMM, file: print_lfmm.c
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
#include "register_lfmm.h"
#include <stdio.h>
#include <stdlib.h>

// print_licence

void print_licence_lfmm()
{
        Rprintf("LFMM Copyright (C) 2012 Eric Frichot\n"
               "This program is free software: you can redistribute it and/or modify\n"
               "it under the terms of the GNU General Public License as published by\n"
               "the Free Software Foundation, either version 3 of the License, or\n"
               "(at your option) any later version.\n"
               "This program is distributed in the hope that it will be useful,\n"
               "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
               "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
               "GNU General Public License for more details.\n"
               "You should have received a copy of the GNU General Public License\n"
               "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n");

}

// print_head_licence

void print_head_licence_lfmm()
{
        Rprintf("LFMM  Copyright (C) 2012 Eric Frichot\n"
               "This program comes with ABSOLUTELY NO WARRANTY; for details type './LFMM -l'.\n"
               "This is free software, and you are welcome to redistribute it\n"
               "under certain conditions; type './LFMM -l' for details.\n\n");

}

// print_head

void print_head_lfmm()
{
        print_head_licence_lfmm();
        Rprintf
            ("****                         LFMM Version 1.3                                   *****\n"
             "****           E. Frichot, S. Schoville, G. Bouchard, O. Francois               *****\n"
             "****                         Please cite our paper !                            *****\n"
             "****   Information at http://membres-timc.imag.fr/Olivier.Francois/lfmm/        *****\n\n");
}

// print_help

void print_help_lfmm()
{

        Rprintf("\nHELP: ./LFMM options \n\n"
               "mandatory:\n"
               "        -x genotype_file        -- genotype file (in .lfmm format)\n"
               "        -v variable_file        -- variable file (in .env format)\n"
               "        -K K                    -- K number of latent factors\n"
               "optional:\n"
               "        -d d                    -- d, the dth covariables         (default: all separately)\n"
               "        -a                      -- all covariables at the same time\n"
               "        -o output_file          -- base of the output files with z-scores (default: genotype_file)\n"
               "        -m                      -- missing data                   (default: no)\n"
               "        -p p                    -- number of processes (CPU)      (default: 1)\n"
               "        -i Niter                -- number of iterations in the GS (default: 1000)\n"
               "        -b burn                 -- burnin parameter in the GS     (default: 500)\n"
               "        -s seed                 -- seed random init               (default: random)\n"
               "        -C dic_file             -- DIC file                       (default: genotype_file.dic)\n"
               "        -h                      -- help\n\n");
}

// print_summary

void print_summary_lfmm(LFMM_param param)
        /*
           int N, int M, int K, int D, int d, int Niter, int burn,
           int m, char *output, char *input, char *cov_file, 
           int num_thrd, long long s, int all)
         */
{

        Rprintf("Summary of the options:\n\n"
               "        -n (number of individuals)      %d\n"
               "        -L (number of loci)             %d\n"
               "        -K (number of latent factors)   %d\n"
               "        -o (output file)                %s\n"
               "        -i (number of iterations)       %d\n"
               "        -b (burnin)                     %d\n"
               "        -s (seed random init)           %llu\n"
               "        -p (number of processes (CPU))  %d\n"
               "        -x (genotype file)              %s\n"
               "        -v (variable file)              %s\n"
               "        -D (number of covariables)      %d\n",
               param->n, param->L, param->K, param->output_file,
               param->Niter, param->burn, param->seed, param->num_thrd,
               param->input_file, param->cov_file, param->D);

        // if missing data
        if (param->nd)
                Rprintf("        -d (the dth covariable)         %d\n",
                       param->nd);
        if (param->all)
                Rprintf("        -a (all variable at the same time)\n");
        if (param->missing_data)
                Rprintf("        -m (missing data)                 \n");

        Rprintf("\n");
}

// print_perc

void print_perc(double *perc, int K, int D)
{
        int k, d;

        Rprintf("\tPercentage of variance:\n");
        for (d = 0; d < D; d++)
                Rprintf("\t\tvar%d\t\t%3.3G %%\n", d, perc[1 + d] * 100);
        for (k = 0; k < K; k++)
                Rprintf("\t\tfactor%d\t\t%3.3G %%\n", k + 1,
                       perc[1 + D + k] * 100);
        Rprintf("\t\tresidual\t%3.3G %%\n", perc[0] * 100);
        Rprintf("\n");

}
