/*
    tracyWidom, file: print_tracyWidom.c
    Copyright (C) 2013 Eric Frichot

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

#include "print_tracyWidom.h"
#include <stdio.h>
#include <stdlib.h>

// print_help

void print_help_tracyWidom()
{
        Rprintf("\nHELP: ./tracyWidom options \n\n"
               "mandatory:\n"
               "        -i input_file         -- eigenvalues file (one column)\n\n"
               "optional:\n"
               "        -h                    -- help\n"
               "        -o output_file        -- output file (default: input_file.tracywidom)\n");
}

// print_summary

void print_summary_tracyWidom(int M, char *input, char *output)
{
        Rprintf("summary of the options:\n\n"
               "        -n (number of eigenvalues)          %d\n"
               "        -i (input file)                     %s\n"
               "        -o (output file)                    %s\n", M, input,
               output);
}
