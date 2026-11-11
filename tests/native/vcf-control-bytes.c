/* Compile the actual private reader with instrumented input/allocation. An
 * invalid prefix has an 8 MiB virtual tail; neither reads nor reallocations
 * may consume it. No large tail is materialized just to test early rejection. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
static const unsigned char *prefix;
static size_t prefix_length, tail_length, position, reads, allocations;
static int probe_fgetc(FILE *file)
{
    (void)file;
    reads++;
    if (position < prefix_length) return prefix[position++];
    if (position - prefix_length < tail_length) { position++; return 'A'; }
    return EOF;
}
static int probe_ferror(FILE *file) { (void)file; return 0; }
static void *probe_realloc(void *ptr, size_t size)
{
    allocations++;
    return realloc(ptr, size);
}
#define fgetc probe_fgetc
#define ferror probe_ferror
#define realloc probe_realloc
#include "convert/vcf2geno.c"
#undef realloc
#undef ferror
#undef fgetc
_Noreturn void error(const char *format, ...)
{
    va_list args; va_start(args, format); vfprintf(stderr, format, args); va_end(args);
    abort();
}
static void input(const unsigned char *bytes, size_t length, size_t tail)
{
    prefix = bytes; prefix_length = length; tail_length = tail;
    position = reads = allocations = 0;
}
int main(void)
{
    int invalid_checks = 0;
    for (int byte = 0; byte <= 127; byte++) {
        if ((byte >= 32 && byte != 127) || byte == '\t' || byte == '\r' || byte == '\n') continue;
        unsigned char first = (unsigned char)byte;
        vcf_context ctx = {0}; ctx.input_name = "invalid-prefix";
        input(&first, 1, 8 * 1024 * 1024);
        assert(vcf_read_line(&ctx) == -1);
        assert(reads == 1 && position == 1 && allocations == 0);
        assert(ctx.line == NULL && ctx.line_capacity == 0 && ctx.message[0]);
        vcf_cleanup(&ctx);
        invalid_checks++;
    }
    /* CR is allowed only as a terminal byte: reject as soon as another
     * non-LF/non-EOF byte follows it, without reading the remaining tail. */
    {
        const unsigned char cr[] = {'\r'};
        vcf_context ctx = {0}; ctx.input_name = "embedded-CR";
        input(cr, sizeof(cr), 8 * 1024 * 1024);
        assert(vcf_read_line(&ctx) == -1);
        assert(reads == 2 && position == 2 && ctx.line_capacity <= 16);
        vcf_cleanup(&ctx);
    }
    const char *valid[] = {"abc\tdef\n", "abc\tdef\r\n", "abc\tdef", "abc\tdef\r"};
    for (size_t i = 0; i < sizeof(valid) / sizeof(valid[0]); i++) {
        vcf_context ctx = {0}; ctx.input_name = "valid-line";
        input((const unsigned char *)valid[i], strlen(valid[i]), 0);
        assert(vcf_read_line(&ctx) == 1 && !strcmp(ctx.line, "abc\tdef"));
        assert(vcf_read_line(&ctx) == 0);
        vcf_cleanup(&ctx);
    }
    printf("%d forbidden-byte prefixes rejected after one read, zero allocations; "
           "embedded CR rejected after two reads; 4 tab/LF/CRLF/EOF cases passed.\n", invalid_checks);
    return 0;
}
