/* Standalone ASan driver; production parsers and I/O are linked unchanged. */
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "convert/ancestrymap.h"

void Rprintf(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
}

_Noreturn void error(const char *format, ...)
{
    if (format) {
        va_list args;
        va_start(args, format);
        vfprintf(stderr, format, args);
        va_end(args);
    }
    fputc('\n', stderr);
    exit(2);
}

int main(int argc, char **argv)
{
    int n = 0, m = 0;
    if (argc != 4) return 3;
    if (!strcmp(argv[1], "line")) {
        char line[2048], name[ANCESTRYMAP_SNP_NAME_CAPACITY];
        int length = atoi(argv[2]), allele = 0, warning = 0;
        if (length < 1 || length > 2000) return 3;
        memset(line, 'a', length);
        strcpy(line + length, " ind 1\n");
        read_line_ancestrymap(line, &allele, name, 1, 1, "native-test", &warning);
        return strlen(name) == (size_t)length && allele == 1 ? 0 : 4;
    }
    if (!strcmp(argv[1], "geno"))
        ancestrymap2geno(argv[2], argv[3], &n, &m);
    else if (!strcmp(argv[1], "lfmm"))
        ancestrymap2lfmm(argv[2], argv[3], &n, &m);
    else return 3;
    return n > 0 && m > 0 ? 0 : 4;
}
