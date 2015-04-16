/*
   LFMM, file: V.c
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
#include "V.h"
#include "data_lfmm.h"
#include "error_lfmm.h"
#include "../matrix/matrix.h"
#include "../matrix/cholesky.h"
#include "../matrix/inverse.h"
#include "../io/read.h"
#include "../matrix/rand.h"

// update_V

void update_V(LFMM_param param, LFMM_GS_param GS_param)
{
        // m_V = U * (I.*(R - C*beta));                         (K,M)
        create_m(param->U, param->dat, param->mC, param->beta,
                 GS_param->m_V, param->L, param->n, param->mD,
                 param->K, param->num_thrd, 0);
#ifdef DEBUG
        print_debug_NaN(GS_param->m_V, param->K, param->L, "m_V");
#endif
        // cov_V = alpha .* eye(K) + alpha_R .* U*U';           (K,K)
        create_inv_cov(GS_param->inv_cov_V, param->alpha_V,
                       param->alpha_R, param->U, param->K, param->n,
                       1);
#ifdef DEBUG
        print_debug_NaN(GS_param->inv_cov_V, param->K, param->K, "inv_cov_V");
#endif
        /*      mu_V = alpha_R .* inv(cov_V) * m_V;                     (K,M)
           for j=1:M
           V(:,j) = mvnrnd(mu_V(:,j),inv(cov_V));
           end                                                             */
        rand_matrix(param->V, GS_param->m_V, GS_param->inv_cov_V,
                    param->alpha_R, param->K, param->L, 1);
#ifdef DEBUG
        print_debug_NaN(param->V, param->K, param->L, "V");
#endif
        // check nan
        if (isnan(param->V[0]))
                print_error_global("nan", NULL, 0);

}
