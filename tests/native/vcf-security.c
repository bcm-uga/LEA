/* Public VCF wrapper driver. Each invocation runs in an isolated process. */
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "R_vcf2geno.h"
void Rprintf(const char *format, ...)
{
    va_list args; va_start(args, format); vfprintf(stdout, format, args); va_end(args);
}
_Noreturn void error(const char *format, ...)
{
    if (format) {
        va_list args; fputs("CONTROLLED_ERROR:", stderr);
        va_start(args, format); vfprintf(stderr, format, args); va_end(args);
    } else fputs("LEGACY_NULL_ERROR", stderr);
    fputc('\n', stderr);
    exit(2);
}
int main(int argc, char **argv)
{
    int n=0, m=0;
    if (argc != 3) return 3;
    char *input=argv[1], *output=argv[2];
    R_vcf2geno(&input, &output, &n, &m);
    printf("VCF_COUNTS:%d:%d\n", n, m);
    return 0;
}
