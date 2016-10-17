/*
    tracyWidom, file: register_tracyWidom.c
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
#include <string.h>
#include <math.h>
#include "register_tracyWidom.h"
#include "error_tracyWidom.h"
#include "print_tracyWidom.h"
#include "../io/io_tools.h"

// analyse_param

void analyse_param_tracyWidom(int argc, char *argv[], char *input, char *output)
{
        int i;
        char *tmp_file;
        int g_data = -1;

        for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
                                // help
                        case 'h':      // global
                                print_help_tracyWidom();
                                error(NULL);
                                break;
                                // input file
                        case 'i':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_tracyWidom("cmd",
                                                               "i (intput file)");
                                g_data = 0;
                                strcpy(input, argv[i]);
                                break;
                                // output file
                        case 'o':
                                i++;
                                if (argc == i || argv[i][0] == '-')
                                        print_error_tracyWidom("cmd",
                                                               "o (output file)");
                                strcpy(output, argv[i]);
                                break;
                        default:
                                print_error_tracyWidom("basic", NULL);
                        }
                } else {
                        print_error_tracyWidom("basic", NULL);
                }
        }

        if (g_data == -1)
                print_error_tracyWidom("option", "-i input_file");

        // write output file names
        tmp_file = remove_ext(input, '.', '/');

        if (!strcmp(output, "")) {
                strcpy(output, tmp_file);
                strcat(output, ".tracywidom");
        }

        Free(tmp_file);
}
