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
#include "als.h"
#include "als_k1.h"
#include "criteria.h"
#include "print_snmf.h"
#include "register_snmf.h"
#include "../matrix/matrix.h"
#include "../matrix/rand.h"
#include "../io/io_data_double.h"
#include "../io/io_tools.h"
#include "../bituint/io_geno_bituint.h"
#include "../bituint/bituint.h"
#include "../createDataSet/createDataSet.h"
#include "../crossEntropy/crossEntropy.h"

void sNMF(sNMF_param param) 
{	
	// temporary variables
	double like = 0.0;
        int K = param->K;
        int n = param->n;
        int L = param->L;
	char *tmp_file; 
	int Mc, Mci, Mp;
	bituint* X;

	//  random init
	init_random(&(param->seed));

	// fix the number of possible factors 
	if (param->m)
		param->nc = param->m + 1;
	else { 
		param->nc = 3;
		param->m = 2;
	}

	// count the number of lines and columns
	param->n = nb_cols_geno(param->input_file);
	param->L = nb_lines(param->input_file, param->n);

	n = param->n;
	L = param->L;
	param->Mc = param->L * param->nc;
	// memory allocation
	param->temp1 = Calloc(K * K * sizeof(double), double);
	param->tempQ = Calloc(n * K * sizeof(double), double);
	param->temp3 = Calloc(n * K * sizeof(double), double);
	param->Y = Calloc(K * n * sizeof(double), double);

	if (param->I == -1) 
		param->I = imin(10000, L/10);
	// write command line summary
        print_summary_snmf(param); 

        // write input file name
	if (param->pourcentage) {
	        tmp_file = remove_ext(param->input_file,'.','/');
                strcpy(param->data_file, tmp_file);
	        strcat(param->data_file, "_I.geno");
	        Free(tmp_file);
		// create file with masked genotypes
		Rprintf("\n <<<<<< createDataSet program\n\n");
		createDataSet(param->input_file, param->seed, param->pourcentage, 
			param->data_file);
		Rprintf("\n >>>>>>\n\n");
	} else 
                strcpy(param->data_file, param->input_file);
	
	// memory allocation
	Mc = param->nc*L;
	init_mat_bituint(&(param->X), n, Mc, &(param->Mp));
        param->Q = (double *) Calloc(n * K * sizeof(double), double);      // of size NxK

	// read of genotypic data
	read_geno_bituint(param->data_file, n, L, param->Mp, param->nc, param->X);
        Rprintf("Read genotype file %s:\t\tOK.\n\n",param->input_file);

	
	// init with a given matrix Q
	if (strcmp(param->input_file_Q,"")) {
		read_data_double(param->input_file_Q, n, K, param->Q);	
	// init of Q with a smaller data set
	} else {
		rand_matrix_double(param->Q, n, K);
		if (param->I && K > 1) {
			// save X
			Mp = param->Mp;
			X = param->X;
			// init subset matrices
			Rprintf("Initialization of Q with a random subset of %d SNPs:\n", param->I);
			Mci = param->nc * param->I;
			init_mat_bituint(&(param->X), n, Mci, &(param->Mp));
        		param->F = (double *) Calloc(K * Mci * sizeof(double), double);     // of size McxK
			// save L
			L = param->L;
			param->L = param->I;
			// save Mc
			Mc = param->Mc;
			param->Mc = Mci;
			// select a subset of SNPs
			select_geno_bituint(X, param->X, n, L, param->I, param->nc, param->Mp, Mp);
			// calc init of Q_res
			ALS(param);
			// free memory
			Free(param->F);
			Free(param->X);
			// put back the parameters
			param->X = X;
			param->Mp = Mp;
			param->L = L;
			param->Mc = Mc; 
		}
	} 

	// memory allocation
        param->F = (double *) Calloc(K * Mc * sizeof(double), double);     // of size McxK

	// parameter estimation
	Rprintf("\nMain algorithm:\n");
	if (K == 1) 
		ALS_k1(param);
	else
		ALS(param);

	// least square estimates
	like = least_square(param); 
	Rprintf("\nLeast-square error: %f\n", like);

	// write Q
	write_data_double(param->output_file_Q, n, K, param->Q);
       	Rprintf("Write individual ancestry coefficient file %s:"
		"\t\tOK.\n",param->output_file_Q);

	// write F
	write_data_double(param->output_file_F, Mc, K ,param->F);
        Rprintf("Write ancestral allele frequency coefficient file %s:"
		"\tOK.\n\n",param->output_file_F);

	// cross-entropy
	if (param->pourcentage) {
		Rprintf("<<<<<< crossEntropy program\n\n");
		crossEntropy(param->input_file, param->data_file, param->output_file_Q, 
		param->output_file_F, K, param->m, &(param->all_ce), &(param->masked_ce));	
		Rprintf("\n >>>>>>\n\n");
	}
}
