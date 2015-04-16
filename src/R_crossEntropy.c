/*
   crossEntropy, file: R_crossEntropy.c
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

#define SEP " "

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "crossEntropy/crossEntropy.h"

#include "R_crossEntropy.h"

void R_crossEntropy(char **R_genotype_file, char **R_missing_data_file,
                    char **R_Q_file, char **R_F_file, int *R_K, int *R_ploidy,
                    double *all_ce, double *masked_ce)
{
        crossEntropy(*R_genotype_file,
                     *R_missing_data_file,
                     *R_Q_file, *R_F_file, *R_K, *R_ploidy, all_ce, masked_ce);
}
