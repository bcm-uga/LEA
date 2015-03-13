/*
    matrix, file: error.c
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

#include "error_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// print_error_global

void print_error_global(char *msg, char *file, int n)
{
        Rprintf("\n");
        // open
        if (!strcmp(msg, "open")) {
                Rprintf
                    ("ERROR: unable to open file %s. Please, check that the name"
                     " of the file you provided is correct.\n", file);
                // read
        } else if (!strcmp(msg, "read")) {
                Rprintf("ERROR: unable to read file %s. Please, check that the"
                       " format is correct (refer to the documentation).\n",
                       file);
                // interne
        } else if (!strcmp(msg, "interne")) {
                Rprintf("ERROR: internal error. Please run the program again. If"
                       " the error is repeated, contact us.\n");
                // constant
        } else if (!strcmp(msg, "constant")) {
                Rprintf("ERROR: %d SNPs are invariant. Please, remove these SNPs"
                       " before the analysis.\n", n);
                // nan
        } else if (!strcmp(msg, "nan")) {
                Rprintf
                    ("ERROR: internal error. Please, run the program again. If"
                     " the error is still present, contact us.\n");
                // else 
        } else {
                Rprintf("ERROR: internal error.\n");
        }

        Rprintf("\n");
        error(NULL);
}
