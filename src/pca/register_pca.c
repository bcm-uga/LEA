/*
    sNMF, file: register_pca.c
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
#include "register_pca.h"
#include "error_pca.h"
#include "print_pca.h"
#include "../io/io_tools.h"

// analyse_param

void analyse_param_pca(	int argc, char *argv[], char *input, char* output_eva, 
			char* output_eve, char* output_sdev, char* output_x,
			int *K, int *c, int *s) 
{
        int i;
	char *tmp_file;
	char tmp[512];
	int g_data = -1;

	for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
			// the number of principal component
                        case 'K':
                                i++;
                                if (argc == i || argv[i][0] == '-')
					print_error_pca("cmd","K (number of principal components)");
                                *K = atoi(argv[i]);
				strcpy(tmp,argv[i]);
                                break;
			// help
                        case 'h':   // global
                                print_help_pca();
                                error(NULL);
                                break;
			// the genotypic file
                        case 'x':
                                i++;
                                if (argc == i || argv[i][0] == '-')
					print_error_pca("cmd","x (genotype file)");
                                g_data = 0;
                                strcpy(input, argv[i]);
                                break;
			// output eigenvalue file
                        case 'a':
                                i++;
                                if (argc == i || argv[i][0] == '-')
					print_error_pca("cmd","a (output eigenvalues file)");
                                strcpy(output_eva, argv[i]);
                                break;
			// output eigenvector file
                        case 'e':
                                i++;
                                if (argc == i || argv[i][0] == '-')
					print_error_pca("cmd","a (output eigenvalues file)");
                                strcpy(output_eve, argv[i]);
                                break;
			// output standard deviation file
                        case 'd':
                                i++;
                                if (argc == i || argv[i][0] == '-')
					print_error_pca("cmd","d (output standard deviation file)");
                                strcpy(output_sdev, argv[i]);
                                break;
			// output projection file
                        case 'p':
                                i++;
                                if (argc == i || argv[i][0] == '-')
					print_error_pca("cmd","p (output projection file)");
                                strcpy(output_x, argv[i]);
                                break;
			// center
                        case 'c':
				*c = 1;
				break;
			// scale
                        case 's':
				*s = 1;
				break;
                        default:    print_error_pca("basic",NULL);
                        }
                } else {
                        print_error_pca("basic",NULL);
		}
        }

	// no input file
        if (g_data == -1)
		print_error_pca("option","-x genotype_file");

	// negative K
        if (*K <= 0)
		*K = 0;

        // write output file names
        tmp_file = remove_ext(input,'.','/');
	if (!strcmp(output_eva,"")) {
		strcpy(output_eva,tmp_file);
		strcat(output_eva,".eigenvalues");
	}
	if (!strcmp(output_eve,"")) {
		strcpy(output_eve,tmp_file);
		strcat(output_eve,".eigenvectors");
	}
	if (!strcmp(output_sdev,"")) {
		strcpy(output_sdev,tmp_file);
		strcat(output_sdev,".sdev");
	}
	if (!strcmp(output_x,"")) {
		strcpy(output_x,tmp_file);
		strcat(output_x,".projection");
	}
	Free(tmp_file);
}

