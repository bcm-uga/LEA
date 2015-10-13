/*
   LFMM, file: beta.c
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
#include "../matrix/matrix.h"
#include "../matrix/cholesky.h"
#include "../matrix/inverse.h"
#include "../matrix/rand.h"
#include "beta.h"
#include "register_lfmm.h"
#include "lfmm_algo.h"
#include "data_lfmm.h"
#include "error_lfmm.h"

// update_beta

void update_beta(LFMM_param param, LFMM_GS_param GS_param)
{
        // m_beta = C' * (R - U'*V);                       (D,M)
        create_m(param->mC, param->dat, param->U, param->V,
                 GS_param->m_beta, param->L, param->n, param->K,
                 param->mD, param->num_thrd, 0);
#ifdef DEBUG
        print_debug_NaN(GS_param->m_beta, param->mD, param->L, "m_beta");
#endif
        // cov_beta = alpha_beta .* eye(D) + alpha_R .* C'*C;   (D,D)
        create_inv_cov(GS_param->inv_cov_beta, param->alpha_beta,
                       param->alpha_R, param->mC, param->mD, param->n,
                       1);
#ifdef DEBUG
        print_debug_NaN(GS_param->inv_cov_beta, param->mD, param->mD,
                        "inv_cov_beta");
#endif
        /*      mu_beta = alpha_R .* inv(cov_beta) * m_beta;            (D,M)
           for j=1:M
           beta(:,j) = mvnrnd(mu_beta(:,j),inv(cov_beta));
           end                                                             */
        rand_matrix(param->beta, GS_param->m_beta, GS_param->inv_cov_beta,
                    param->alpha_R, param->mD, param->L, 1);
#ifdef DEBUG
        print_debug_NaN(param->beta, param->mD, param->L, "beta");
#endif
        // nan check
        if (isnan(param->beta[0]))
                print_error_global("nan", NULL, 0);
}

// update_alpha_beta

void update_alpha_beta(LFMM_param param)
{
        // temporary parameters
        double *beta = param->beta;
        double b_epsilon = param->b_epsilon;
        int D = param->mD;
        int M = param->L;
        int d;

        //int a = (int)(epsilon) + M/2;
        int a = (int)1 + M / 2;
        // allocate memory
        double *bb = (double *)Calloc(D * sizeof(double), double);

        // b = 1/2 * sum(sum(beta.^2));
        dble_sum2(beta, D, M, bb, b_epsilon);   // beta(K,M)
        //dble_sum2(beta,D,M,bb,1); // beta(K,M)

        // update alpha_beta
        //param->alpha_beta[0] = rand_gamma(a+ b_epsilon, 1.0/(double)(bb[0]+ b_epsilon));
        param->alpha_beta[0] =
            rand_gamma(a + b_epsilon - 1.0,
                       1.0 / (double)(bb[0] + b_epsilon - 1.0));
        for (d = 1; d < D; d++)
                param->alpha_beta[d] = rand_gamma(a, 1.0 / (double)(bb[d]));

        // free memory
        Free(bb);
}
