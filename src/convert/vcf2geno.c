/*
 *     vcf2geno, file: vcf2geno.c
 *     Copyright (C) 2013 Eric Frichot
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <limits.h>
#include "vcf2geno.h"

/* All storage and streams belong to this context. Helpers return errors rather
 * than jumping into R, so cleanup always precedes the R error. Field pointers
 * borrow the complete line and are used only until the next read. */
typedef struct {
        const char *input_name;
        FILE *input, *geno, *snp, *removed;
        char *line;
        size_t line_capacity, line_number;
        char **fields, **keys;
        size_t field_capacity, field_count, key_capacity;
        int *allele;
        char message[1024];
} vcf_context;

static int vcf_fail(vcf_context *ctx, const char *format, ...)
{
        char detail[512];
        va_list args;
        va_start(args, format);
        vsnprintf(detail, sizeof(detail), format, args);
        va_end(args);
        snprintf(ctx->message, sizeof(ctx->message),
                 "VCF error at line %zu: %s (file: %.400s)",
                 ctx->line_number, detail, ctx->input_name);
        return 0;
}

static void *vcf_grow(vcf_context *ctx, void *old, size_t *capacity,
                      size_t needed, size_t width)
{
        size_t next;
        void *replacement;
        if (needed <= *capacity)
                return old;
        if (needed > SIZE_MAX / width) {
                vcf_fail(ctx, "Buffer size exceeds allocation limits.");
                return NULL;
        }
        next = *capacity ? *capacity : 16;
        while (next < needed) {
                if (next > SIZE_MAX / 2) {
                        next = needed;
                        break;
                }
                next *= 2;
        }
        if (next > SIZE_MAX / width)
                next = needed;
        replacement = realloc(old, next * width);
        if (!replacement) {
                vcf_fail(ctx, "Unable to allocate a VCF buffer.");
                return NULL;
        }
        *capacity = next;
        return replacement;
}

/* 1 = a complete physical line, 0 = clean EOF, -1 = error. Use int for EOF,
 * binary input plus explicit CRLF handling, and no POSIX-only getline API. */
static int vcf_read_line(vcf_context *ctx)
{
        size_t length = 0;
        int byte = fgetc(ctx->input);
        if (byte == EOF) {
                if (ferror(ctx->input)) {
                        vcf_fail(ctx, "Unable to read input.");
                        return -1;
                }
                return 0;
        }
        if (ctx->line_number == SIZE_MAX) {
                vcf_fail(ctx, "Too many physical lines.");
                return -1;
        }
        ctx->line_number++;
        while (byte != EOF && byte != '\n') {
                char *replacement;
                if (byte == 0) {
                        vcf_fail(ctx, "Embedded NUL in input.");
                        return -1;
                }
                if ((byte < 32 && byte != '\t' && byte != '\r') || byte == 127 ||
                    (length && ctx->line[length - 1] == '\r')) {
                        vcf_fail(ctx, "Invalid control byte in input.");
                        return -1;
                }
                if (length > SIZE_MAX - 2) {
                        vcf_fail(ctx, "Physical line exceeds allocation limits.");
                        return -1;
                }
                replacement = vcf_grow(ctx, ctx->line, &ctx->line_capacity,
                                       length + 2, sizeof(char));
                if (!replacement)
                        return -1;
                ctx->line = replacement;
                ctx->line[length++] = (char)byte;
                byte = fgetc(ctx->input);
        }
        if (ferror(ctx->input)) {
                vcf_fail(ctx, "Unable to read input.");
                return -1;
        }
        if (!ctx->line) {
                ctx->line = vcf_grow(ctx, NULL, &ctx->line_capacity, 1, sizeof(char));
                if (!ctx->line)
                        return -1;
        }
        if (length && ctx->line[length - 1] == '\r')
                length--;
        ctx->line[length] = '\0';
        return 1;
}

/* Split spaces within each tab-delimited segment for historical LEA inputs,
 * but never collapse an empty tab segment into the next field. */
static int vcf_split(vcf_context *ctx)
{
        char *segment = ctx->line;
        int has_tabs = strchr(segment, '\t') != NULL;
        ctx->field_count = 0;
        for (;;) {
                char *end = strchr(segment, '\t');
                char *cursor = segment;
                size_t before = ctx->field_count;
                if (end)
                        *end = '\0';
                while (*cursor) {
                        char *field;
                        char **replacement;
                        while (*cursor == ' ')
                                cursor++;
                        if (!*cursor)
                                break;
                        field = cursor;
                        while (*cursor && *cursor != ' ')
                                cursor++;
                        if (*cursor)
                                *cursor++ = '\0';
                        if (ctx->field_count == SIZE_MAX)
                                return vcf_fail(ctx, "Too many fields.");
                        replacement = vcf_grow(ctx, ctx->fields, &ctx->field_capacity,
                                               ctx->field_count + 1, sizeof(char *));
                        if (!replacement)
                                return 0;
                        ctx->fields = replacement;
                        ctx->fields[ctx->field_count++] = field;
                }
                if (has_tabs && before == ctx->field_count)
                        return vcf_fail(ctx, "Empty tab-separated field.");
                if (!end)
                        break;
                segment = end + 1;
        }
        return 1;
}

static int compare_strings(const void *left, const void *right)
{
        return strcmp(*(char *const *)left, *(char *const *)right);
}

static int vcf_read_header(vcf_context *ctx, int *samples)
{
        static const char *const names[9] = {
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
        };
        int status;
        while ((status = vcf_read_line(ctx)) == 1) {
                if (!strncmp(ctx->line, "##", 2))
                        continue;
                if (!vcf_split(ctx))
                        return 0;
                if (ctx->field_count < 10)
                        return vcf_fail(ctx, "Expected #CHROM header with FORMAT and at least one sample.");
                for (size_t i = 0; i < 9; i++)
                        if (strcmp(ctx->fields[i], names[i]))
                                return vcf_fail(ctx, "Invalid VCF header column %zu; expected %s.", i + 1, names[i]);
                if (ctx->field_count - 9 > INT_MAX)
                        return vcf_fail(ctx, "Sample count exceeds integer limits.");
                /* Header sample pointers are not needed after this count. */
                qsort(ctx->fields + 9, ctx->field_count - 9, sizeof(char *), compare_strings);
                for (size_t i = 10; i < ctx->field_count; i++)
                        if (!strcmp(ctx->fields[i - 1], ctx->fields[i]))
                                return vcf_fail(ctx, "Duplicate sample identifier in header.");
                *samples = (int)(ctx->field_count - 9);
                return 1;
        }
        if (status < 0)
                return 0;
        return vcf_fail(ctx, "Missing #CHROM header (empty or incomplete VCF input).");
}

static int key_initial(unsigned char c)
{
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
}

/* Temporarily split FORMAT for grammar/duplicate validation, then restore its
 * original bytes before sidecar output. key_count also bounds sample suffixes. */
static int vcf_format(vcf_context *ctx, size_t *key_count)
{
        char *format = ctx->fields[8], *cursor = format;
        size_t length = strlen(format);
        *key_count = 0;
        if (strcspn(format, ":") != 2 || format[0] != 'G' || format[1] != 'T')
                return vcf_fail(ctx, "Unsupported FORMAT: GT must be the first field.");
        for (;;) {
                char *end = strchr(cursor, ':');
                char **replacement;
                if (end)
                        *end = '\0';
                if (!key_initial((unsigned char)*cursor))
                        return vcf_fail(ctx, "Empty or malformed FORMAT key.");
                for (char *p = cursor + 1; *p; p++)
                        if (!key_initial((unsigned char)*p) && !(*p >= '0' && *p <= '9') && *p != '.')
                                return vcf_fail(ctx, "Malformed FORMAT key.");
                if (*key_count == SIZE_MAX)
                        return vcf_fail(ctx, "Too many FORMAT keys.");
                replacement = vcf_grow(ctx, ctx->keys, &ctx->key_capacity,
                                       *key_count + 1, sizeof(char *));
                if (!replacement)
                        return 0;
                ctx->keys = replacement;
                ctx->keys[(*key_count)++] = cursor;
                if (!end)
                        break;
                cursor = end + 1;
        }
        qsort(ctx->keys, *key_count, sizeof(char *), compare_strings);
        for (size_t i = 1; i < *key_count; i++)
                if (!strcmp(ctx->keys[i - 1], ctx->keys[i]))
                        return vcf_fail(ctx, "Duplicate FORMAT key.");
        for (size_t i = 0; i < length; i++)
                if (!format[i])
                        format[i] = ':';
        return 1;
}

static int vcf_allele_char(char c)
{
        return c == '0' || c == '1' || c == '.';
}

static int vcf_decode(vcf_context *ctx, int samples, size_t key_count)
{
        for (size_t i = 0; i < (size_t)samples; i++) {
                const char *sample;
                size_t length, values = 1;
                int allele;
                if (i + 9 >= ctx->field_count)
                        return vcf_fail(ctx, "Missing sample %zu.", i + 1);
                sample = ctx->fields[i + 9];
                length = strcspn(sample, ":");
                if (length == 1 && vcf_allele_char(sample[0]))
                        allele = sample[0] == '.' ? 9 : sample[0] - '0';
                else if (length == 3 && vcf_allele_char(sample[0]) &&
                         (sample[1] == '/' || sample[1] == '|') && vcf_allele_char(sample[2]))
                        allele = (sample[0] == '.' || sample[2] == '.') ? 9 :
                                 (sample[0] - '0') + (sample[2] - '0');
                else
                        return vcf_fail(ctx, "Unsupported or malformed haploid/diploid GT at sample %zu.", i + 1);
                sample += length;
                while (*sample) {
                        sample++; /* colon */
                        length = strcspn(sample, ":");
                        if (!length)
                                return vcf_fail(ctx, "Empty FORMAT value at sample %zu.", i + 1);
                        if (values == key_count)
                                return vcf_fail(ctx, "Too many FORMAT values at sample %zu.", i + 1);
                        values++;
                        sample += length;
                }
                ctx->allele[i] = allele;
        }
        return 1;
}

static int vcf_write_info(vcf_context *ctx, FILE *file, int removed)
{
        for (size_t i = 0; i < 9; i++)
                if (fprintf(file, "%s ", ctx->fields[i]) < 0)
                        return vcf_fail(ctx, "Unable to write variant information.");
        if ((removed ? fputs("REMOVED\n", file) : fputc('\n', file)) == EOF)
                return vcf_fail(ctx, "Unable to write variant information.");
        return 1;
}

static void vcf_close(vcf_context *ctx, FILE **stream, const char *role)
{
        if (*stream) {
                int result = fclose(*stream);
                *stream = NULL;
                if (result == EOF && !ctx->message[0])
                        vcf_fail(ctx, "Unable to close %s.", role);
        }
}

static void vcf_cleanup(vcf_context *ctx)
{
        vcf_close(ctx, &ctx->input, "input");
        vcf_close(ctx, &ctx->geno, "GENO output");
        vcf_close(ctx, &ctx->snp, "SNP information output");
        vcf_close(ctx, &ctx->removed, "removed-variant output");
        free(ctx->allele);
        free(ctx->keys);
        free(ctx->fields);
        free(ctx->line);
}

int nb_cols_vcf(char *file)
{
        vcf_context ctx = {0};
        int samples = 0;
        ctx.input_name = file ? file : "<null>";
        ctx.input = file ? fopen(file, "rb") : NULL;
        if (!ctx.input)
                vcf_fail(&ctx, "Unable to open VCF input.");
        else
                vcf_read_header(&ctx, &samples);
        vcf_cleanup(&ctx);
        if (ctx.message[0])
                error("%s", ctx.message);
        return samples;
}

/* The wrapper owns its paths. Return diagnostics only after releasing every
 * parser resource so the wrapper can free its own allocations before R errors. */
int vcf2geno_checked(char *input_file, char *output_file, int *N, int *M,
                     char *snp_file, char *removed_file, int *removed,
                     char *message, size_t message_size)
{
        vcf_context ctx = {0};
        int samples = 0, kept = 0, excluded = 0, status;
        ctx.input_name = input_file ? input_file : "<null>";
        if (!input_file || !output_file || !snp_file || !removed_file || !N || !M || !removed) {
                vcf_fail(&ctx, "Missing conversion argument.");
                goto cleanup;
        }
        ctx.input = fopen(input_file, "rb");
        if (!ctx.input) {
                vcf_fail(&ctx, "Unable to open VCF input.");
                goto cleanup;
        }
        if (!vcf_read_header(&ctx, &samples))
                goto cleanup;
        if ((size_t)samples > SIZE_MAX / sizeof(int)) {
                vcf_fail(&ctx, "Sample array exceeds allocation limits.");
                goto cleanup;
        }
        ctx.allele = calloc((size_t)samples, sizeof(int));
        if (!ctx.allele) {
                vcf_fail(&ctx, "Unable to allocate genotype array.");
                goto cleanup;
        }
        ctx.snp = fopen(snp_file, "w");
        if (!ctx.snp) {
                vcf_fail(&ctx, "Unable to open SNP information output: %.300s", snp_file);
                goto cleanup;
        }
        ctx.removed = fopen(removed_file, "w");
        if (!ctx.removed) {
                vcf_fail(&ctx, "Unable to open removed-variant output: %.300s", removed_file);
                goto cleanup;
        }
        ctx.geno = fopen(output_file, "w");
        if (!ctx.geno) {
                vcf_fail(&ctx, "Unable to open GENO output: %.300s", output_file);
                goto cleanup;
        }
        while ((status = vcf_read_line(&ctx)) == 1) {
                int is_removed;
                size_t key_count;
                if (ctx.line[0] == '#') {
                        vcf_fail(&ctx, "Unexpected header after #CHROM.");
                        goto cleanup;
                }
                if (!vcf_split(&ctx))
                        goto cleanup;
                if (ctx.field_count < 9 || ctx.field_count - 9 != (size_t)samples) {
                        vcf_fail(&ctx, "Expected nine fixed fields and %d samples; observed %zu fields.",
                                 samples, ctx.field_count);
                        goto cleanup;
                }
                is_removed = strlen(ctx.fields[3]) > 1 || strlen(ctx.fields[4]) > 1;
                if (is_removed) {
                        if (excluded == INT_MAX) {
                                vcf_fail(&ctx, "Removed-variant count exceeds integer limits.");
                                goto cleanup;
                        }
                        if (!vcf_write_info(&ctx, ctx.removed, 1))
                                goto cleanup;
                        excluded++;
                } else {
                        if (kept == INT_MAX) {
                                vcf_fail(&ctx, "Locus count exceeds integer limits.");
                                goto cleanup;
                        }
                        if (!vcf_format(&ctx, &key_count) || !vcf_decode(&ctx, samples, key_count))
                                goto cleanup;
                        if (!vcf_write_info(&ctx, ctx.snp, 0))
                                goto cleanup;
                        for (int i = 0; i < samples; i++)
                                if (fputc('0' + ctx.allele[i], ctx.geno) == EOF) {
                                        vcf_fail(&ctx, "Unable to write GENO output.");
                                        goto cleanup;
                                }
                        if (fputc('\n', ctx.geno) == EOF) {
                                vcf_fail(&ctx, "Unable to write GENO output.");
                                goto cleanup;
                        }
                        kept++;
                }
        }
        if (status < 0)
                goto cleanup;
cleanup:
        vcf_cleanup(&ctx);
        if (message && message_size)
                snprintf(message, message_size, "%s", ctx.message);
        if (ctx.message[0])
                return 0;
        *N = samples;
        *M = kept;
        *removed = excluded;
        return 1;
}

void vcf2geno(char *input_file, char *output_file, int *N, int *M,
              char *snp_file, char *removed_file, int *removed)
{
        char message[1024];
        if (!vcf2geno_checked(input_file, output_file, N, M, snp_file,
                              removed_file, removed, message, sizeof(message)))
                error("%s", message);
}
