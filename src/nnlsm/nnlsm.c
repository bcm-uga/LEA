/*
   NMF, file: nnlsm.c
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nnlsm.h"
#include "blockpivot.h"
#include "solvenormaleqcomb.h"
#include "../matrix/inverse.h"

// allocate_nnlsm

Nnlsm_param allocate_nnlsm(int N, int K)
{
	Nnlsm_param param = (Nnlsm_param) Calloc(1 *  sizeof(nnlsm_param), nnlsm_param);

	param->P = (int *)Calloc(N * sizeof(int), int);
	param->Ninf = (int *)Calloc(N * sizeof(int), int);
	param->PassiveSet = (int *)Calloc(N*K * sizeof(int), int);
	param->NonOptSet = (int *)Calloc(N*K * sizeof(int), int);
	param->InfeaSet = (int *)Calloc(N*K * sizeof(int), int);
	param->NotGood = (int *)Calloc(N * sizeof(int), int);
	param->Cols3Ix = (int *)Calloc(N * sizeof(int), int);
	param->subX = (double *)Calloc(N*K * sizeof(double), double);
	param->subY = (double *)Calloc(N*K * sizeof(double), double);
	param->subAtB = (double *)Calloc(N*K * sizeof(double), double);
	param->subPassiveSet = (int *)Calloc(N*K * sizeof(int), int);
	param->selectK = (int *)Calloc(K * sizeof(int), int);
	param->selectN = (int *)Calloc(N * sizeof(int), int);
	param->breaks = (int *)Calloc(N * sizeof(int), int);
	param->sortIx = (int *)Calloc(N * sizeof(int), int);
	param->sAtA = (double *)Calloc(K*K * sizeof(double), double);
	param->inVsAtA = (double *)Calloc(K*K * sizeof(double), double);
        param->tempSortIx = (int *)Calloc(N * sizeof(int), int);
        param->Y = (double *)Calloc(K*N * sizeof(double), double);

	return param;
}

// free_nnlsm

void free_nnlsm(Nnlsm_param param)
{
	Free(param->P);
	Free(param->Ninf);
	Free(param->PassiveSet);
	Free(param->NonOptSet);
	Free(param->InfeaSet);
	Free(param->NotGood);
	Free(param->Cols3Ix);
	Free(param->subX);
	Free(param->subY);
	Free(param->subAtB);
	Free(param->subPassiveSet);
	Free(param->selectK);
	Free(param->selectN);
	Free(param->breaks);
	Free(param->sortIx);
	Free(param->sAtA);
	Free(param->inVsAtA);
	Free(param->tempSortIx);
	Free(param->Y);
}

// nnlsm_blockpivot

int nnlsm_blockpivot(double* AtA, double* AtB, int N, int K, double *X, double *Y,
	Nnlsm_param param)
{
	int niter = 0, bigiter, zeros;
	int maxiter = 5*K;
	int* P = param->P;
	int* Ninf = param->Ninf;
	int* PassiveSet = param->PassiveSet;
	int* NonOptSet = param->NonOptSet;
	int* InfeaSet = param->InfeaSet;
	int* NotGood = param->NotGood;
	int* Cols3Ix = param->Cols3Ix;
	int NotOptCols_nb, pbar;

	// init
	pbar = 3;
	zeros = parameter_init(PassiveSet, NotGood, Ninf, P, K, N, X);
	if (!zeros)
		niter += XY_update(AtA, AtB, PassiveSet, NotGood, N, N, K, X, Y, param);
	else 
		update_Y(AtA, AtB, X, Y, N, K);
	// opt_param
	opt_param_update(PassiveSet, NotGood, NonOptSet, InfeaSet, X, Y,&NotOptCols_nb, N, K);
	// main loop
	bigiter = 1;
	while (NotOptCols_nb && bigiter <= maxiter) {
		bigiter ++;
		// P and PassiveSet update
		PassiveSet_update(NotGood, Ninf, P, pbar, NonOptSet, PassiveSet, InfeaSet, N, K, Cols3Ix);
		// X and Y update
		niter += XY_update(AtA, AtB, PassiveSet, NotGood, NotOptCols_nb, N, K, X, Y, param);
		// opt_param update
		opt_param_update(PassiveSet, NotGood, NonOptSet, InfeaSet, X , Y, &NotOptCols_nb, N, K);
		bigiter++;
	}

	return niter;
}

