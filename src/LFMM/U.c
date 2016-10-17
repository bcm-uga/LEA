/*
   LFMM, file: U.c
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
#include "U.h"
#include "data_lfmm.h"
#include "../matrix/matrix.h"
#include "../matrix/cholesky.h"
#include "../matrix/inverse.h"
#include "../io/read.h"
#include "../matrix/rand.h"
#include "error_lfmm.h"

// update_U

void update_U(LFMM_param param, LFMM_GS_param GS_param)
{
        // m_U = (I.*(R - C*beta)) * V';                        (N,K)
        create_m(param->V, param->dat, param->mC, param->beta,
                 GS_param->m_U, param->L, param->n, param->mD,
                 param->K, param->num_thrd, 1);
#ifdef DEBUG
        print_debug_NaN(GS_param->m_U, param->K, param->n, "m_U");
#endif
        // cov_U = alpha .* eye(K) + alpha_R .* V*V';           (K,K)
        create_inv_cov(GS_param->inv_cov_U, param->alpha_U, param->alpha_R,
                       param->V, param->K, param->L, param->num_thrd);
#ifdef DEBUG
        print_debug_NaN(GS_param->inv_cov_U, param->K, param->K, "inv_cov_U");
#endif
        /*      mu_U = alpha_R .* inv(cov_U) * m_U';                    (N,K)
           for i=1:N
           U(:,i) = mvnrnd(mu_U(:,i),inv(cov_U));
           end                                                             */
        rand_matrix(param->U, GS_param->m_U, GS_param->inv_cov_U,
                    param->alpha_R, param->K, param->n, 1);
#ifdef DEBUG
        print_debug_NaN(param->U, param->K, param->n, "U");
#endif
        if (isnan(param->U[0]))
                print_error_global("nan", NULL, 0);
}

// update_alpha_U

void update_alpha_U(LFMM_param param)
{
        // temporary parameters
        int N = param->n;
        int K = param->K;
        double epsilon = param->b_epsilon;
        double *U = param->U;
        int k;

        // int a = (int)epsilon + N * K / 2;
        int a = (int)epsilon + N / 2;
        // b = 1/2*sum(sum(U.^2)) + 1/2*sum(sum(V.^2));
        // WARNING: ceci est une erreur dans le code (ca devrait etre un double et non un int)
        int b;                  // = 0.5 * dble_sum(U, N * K) + epsilon; // U(D,N) V(D,M)

        // update alpha_U
        //tmp = rand_gamma(a, 1.0 / b);
        for (k = 0; k < K; k++) {
                b = 0.5 * dble_sum(&(U[k * N]), N) + epsilon;
                param->alpha_U[k] = rand_gamma(a, 1.0 / b);
                param->alpha_V[k] = 1.0;
        }
}
