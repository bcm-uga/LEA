/*
 *     convert, file: ancestrymap.c
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
#include <errno.h>
#include "ancestrymap.h"

/* Own streams/storage here; checked helpers never jump into R. The legacy
 * line and matrix reader entry points may borrow caller storage explicitly. */
typedef struct {
        const char *input_name;
        FILE *input, *output;
        char *line;
        size_t capacity, line_number, cells;
        char reference[ANCESTRYMAP_SNP_NAME_CAPACITY];
        int *data, owns_data, N, M, warned;
        size_t warning_i, warning_j;
        char message[1024];
} am_context;

static int am_fail(am_context *ctx, const char *format, ...)
{
        char detail[512];
        va_list args;
        va_start(args, format);
        vsnprintf(detail, sizeof(detail), format, args);
        va_end(args);
        snprintf(ctx->message, sizeof(ctx->message), "%s (file: %.350s, line %zu)",
                 detail, ctx->input_name, ctx->line_number);
        return 0;
}

static void am_close(am_context *ctx, FILE **stream, const char *role)
{
        if (*stream) {
                int status = fclose(*stream);
                *stream = NULL;
                if (status == EOF && !ctx->message[0])
                        am_fail(ctx, "Unable to close ANCESTRYMAP %s.", role);
        }
}

static void am_cleanup(am_context *ctx)
{
        am_close(ctx, &ctx->input, "input");
        am_close(ctx, &ctx->output, "output");
        if (ctx->owns_data)
                free(ctx->data);
        free(ctx->line);
}

static void am_warning(const am_context *ctx)
{
        if (ctx->warned) {
                Rprintf("Warning: some genotypes are not 0, 1, 2 or 9 in %s.\n", ctx->input_name);
                /* Preserve the historical warning coordinate convention. */
                Rprintf("\t First warning at individual %zu, column %zu.\n\n",
                        ctx->warning_j, ctx->warning_i);
        }
}

static int am_noncanonical(int value)
{
        return value != 0 && value != 1 && value != 2 && value != 9;
}

static int am_reserve(am_context *ctx, size_t needed)
{
        size_t capacity = ctx->capacity ? ctx->capacity : 256;
        char *replacement;
        if (needed <= ctx->capacity)
                return 1;
        while (capacity < needed) {
                if (capacity > SIZE_MAX / 2) {
                        capacity = needed;
                        break;
                }
                capacity *= 2;
        }
        replacement = realloc(ctx->line, capacity);
        if (!replacement)
                return am_fail(ctx, "Unable to allocate ANCESTRYMAP line buffer.");
        ctx->line = replacement;
        ctx->capacity = capacity;
        return 1;
}

/* Complete physical lines, including final lines without LF. 1=line, 0=EOF,
 * -1=error. Standard C only; binary input makes CRLF handling explicit. */
static int am_read_line(am_context *ctx)
{
        size_t length = 0;
        int byte = fgetc(ctx->input);
        if (byte == EOF) {
                if (ferror(ctx->input)) {
                        am_fail(ctx, "Error reading ANCESTRYMAP input.");
                        return -1;
                }
                return 0;
        }
        if (ctx->line_number == SIZE_MAX) {
                am_fail(ctx, "ANCESTRYMAP physical line count exceeds limits.");
                return -1;
        }
        ctx->line_number++;
        while (byte != EOF && byte != '\n') {
                if (byte == 0 || (byte < 32 && byte != '\t' && byte != '\r') || byte == 127 ||
                    (length && ctx->line[length - 1] == '\r')) {
                        am_fail(ctx, "Invalid ANCESTRYMAP record: embedded NUL or control byte.");
                        return -1;
                }
                if (length > SIZE_MAX - 2) {
                        am_fail(ctx, "ANCESTRYMAP line exceeds allocation limits.");
                        return -1;
                }
                if (!am_reserve(ctx, length + 2))
                        return -1;
                ctx->line[length++] = (char)byte;
                byte = fgetc(ctx->input);
        }
        if (ferror(ctx->input)) {
                am_fail(ctx, "Error reading ANCESTRYMAP input.");
                return -1;
        }
        if (!am_reserve(ctx, length + 1))
                return -1;
        if (length && ctx->line[length - 1] == '\r')
                length--;
        ctx->line[length] = '\0';
        return 1;
}

static int am_space(char byte)
{
        return byte == ' ' || byte == '\t' || byte == '\r' || byte == '\n';
}

static char *am_token(char **cursor)
{
        char *token;
        while (**cursor && am_space(**cursor))
                (*cursor)++;
        if (!**cursor)
                return NULL;
        token = *cursor;
        while (**cursor && !am_space(**cursor))
                (*cursor)++;
        if (**cursor)
                *(*cursor)++ = '\0';
        return token;
}

/* Ignore historical annotation fields after the third token, but only after
 * the physical reader has validated the entire line as non-NUL text. */
static int am_parse(am_context *ctx, char *text, char **name, int *allele)
{
        char *cursor = text, *token, *end;
        long value;
        if (!text || !(*name = am_token(&cursor)))
                return am_fail(ctx, "Invalid ANCESTRYMAP record: missing SNP name.");
        if (strlen(*name) >= ANCESTRYMAP_SNP_NAME_CAPACITY)
                return am_fail(ctx, "SNP name too long: maximum is %d bytes.", ANCESTRYMAP_SNP_NAME_CAPACITY - 1);
        if (!am_token(&cursor))
                return am_fail(ctx, "Invalid ANCESTRYMAP record: missing individual identifier and genotype.");
        token = am_token(&cursor);
        if (!token)
                return am_fail(ctx, "Invalid ANCESTRYMAP record: missing genotype.");
        errno = 0;
        value = strtol(token, &end, 10);
        if (end == token || *end || errno == ERANGE || value < INT_MIN || value > INT_MAX)
                return am_fail(ctx, "Invalid ANCESTRYMAP record: invalid or out-of-range integer genotype.");
        *allele = (int)value;
        return 1;
}

static int am_increment(am_context *ctx, size_t *records)
{
        if (*records >= (size_t)INT_MAX)
                return am_fail(ctx, "ANCESTRYMAP record count exceeds integer limits.");
        (*records)++;
        return 1;
}

static int am_count(am_context *ctx, int first_only)
{
        size_t records = 0, individuals = 0;
        int status, changed = 0, allele;
        char *name;
        while ((status = am_read_line(ctx)) == 1) {
                if (!am_parse(ctx, ctx->line, &name, &allele))
                        return 0;
                if (!records)
                        memcpy(ctx->reference, name, strlen(name) + 1);
                if (!changed && strcmp(ctx->reference, name)) {
                        changed = 1;
                        if (first_only)
                                break;
                }
                if (!am_increment(ctx, &records))
                        return 0;
                if (!changed)
                        individuals++;
        }
        if (status < 0)
                return 0;
        if (!records)
                return am_fail(ctx, "Empty ANCESTRYMAP file: no records.");
        ctx->N = (int)individuals; /* bounded by the checked record count */
        if (!first_only) {
                if (records % individuals)
                        return am_fail(ctx, "Invalid ANCESTRYMAP file: incorrect number of lines for %d individuals.", ctx->N);
                ctx->M = (int)(records / individuals);
                ctx->cells = records;
        }
        return 1;
}

static int am_matrix_size(am_context *ctx)
{
        if (ctx->N <= 0 || ctx->M <= 0 || (size_t)ctx->N > SIZE_MAX / (size_t)ctx->M)
                return am_fail(ctx, "Invalid ANCESTRYMAP matrix dimensions.");
        ctx->cells = (size_t)ctx->N * (size_t)ctx->M;
        if (ctx->cells > (size_t)INT_MAX || ctx->cells > SIZE_MAX / sizeof(int))
                return am_fail(ctx, "ANCESTRYMAP matrix exceeds integer or allocation limits.");
        return 1;
}

static int am_read_matrix(am_context *ctx)
{
        size_t i = 0, j = 0, records = 0;
        int status, allele;
        char *name;
        while ((status = am_read_line(ctx)) == 1) {
                if (!am_parse(ctx, ctx->line, &name, &allele))
                        return 0;
                if (!ctx->warned && am_noncanonical(allele)) {
                        ctx->warned = 1;
                        ctx->warning_i = i + 1;
                        ctx->warning_j = j + 1;
                }
                if (records >= ctx->cells)
                        return am_fail(ctx, "Invalid ANCESTRYMAP group sizes: unexpected extra records.");
                if (!records)
                        memcpy(ctx->reference, name, strlen(name) + 1);
                if (strcmp(ctx->reference, name)) {
                        if (i != (size_t)ctx->N)
                                return am_fail(ctx, "Invalid ANCESTRYMAP file at SNP %zu: %zu individuals instead of %d.", j + 1, i, ctx->N);
                        i = 0;
                        j++;
                        memcpy(ctx->reference, name, strlen(name) + 1);
                }
                if (i >= (size_t)ctx->N || j >= (size_t)ctx->M)
                        return am_fail(ctx, "Invalid ANCESTRYMAP group sizes at individual %zu, SNP %zu.", i + 1, j + 1);
                ctx->data[i * (size_t)ctx->M + j] = allele;
                i++;
                records++;
        }
        if (status < 0)
                return 0;
        if (records != ctx->cells || i != (size_t)ctx->N || j + 1 != (size_t)ctx->M)
                return am_fail(ctx, "Invalid ANCESTRYMAP file: inconsistent number of individuals or loci.");
        return 1;
}

static int am_write(am_context *ctx, int lfmm)
{
        size_t rows = (size_t)(lfmm ? ctx->N : ctx->M);
        size_t columns = (size_t)(lfmm ? ctx->M : ctx->N);
        for (size_t row = 0; row < rows; row++) {
                for (size_t column = 0; column < columns; column++) {
                        size_t index = lfmm ? row * (size_t)ctx->M + column : column * (size_t)ctx->M + row;
                        if (fprintf(ctx->output, "%d", ctx->data[index]) < 0 ||
                            (lfmm && column + 1 < columns && fputc(' ', ctx->output) == EOF))
                                return am_fail(ctx, "Unable to write ANCESTRYMAP conversion output.");
                }
                if (fputc('\n', ctx->output) == EOF)
                        return am_fail(ctx, "Unable to write ANCESTRYMAP conversion output.");
        }
        return 1;
}

static void am_convert(char *input, char *output, int *N, int *M, int lfmm)
{
        am_context ctx = {0};
        ctx.input_name = input ? input : "<null>";
        if (!input || !*input || !output || !*output || !N || !M) {
                am_fail(&ctx, "Invalid ANCESTRYMAP conversion arguments.");
                goto cleanup;
        }
        ctx.input = fopen(input, "rb");
        if (!ctx.input) {
                am_fail(&ctx, "Unable to open ANCESTRYMAP input.");
                goto cleanup;
        }
        if (!am_count(&ctx, 0) || !am_matrix_size(&ctx))
                goto cleanup;
        ctx.data = calloc(ctx.cells, sizeof(int));
        ctx.owns_data = 1;
        if (!ctx.data) {
                am_fail(&ctx, "Unable to allocate ANCESTRYMAP matrix.");
                goto cleanup;
        }
        if (fseek(ctx.input, 0, SEEK_SET)) {
                am_fail(&ctx, "Unable to rewind ANCESTRYMAP input.");
                goto cleanup;
        }
        ctx.line_number = 0;
        if (!am_read_matrix(&ctx))
                goto cleanup;
        am_close(&ctx, &ctx.input, "input");
        if (ctx.message[0])
                goto cleanup;
        ctx.output = fopen(output, "w");
        if (!ctx.output) {
                am_fail(&ctx, "Unable to open ANCESTRYMAP output: %.250s", output);
                goto cleanup;
        }
        if (!am_write(&ctx, lfmm))
                goto cleanup;
cleanup:
        am_cleanup(&ctx);
        am_warning(&ctx);
        if (ctx.message[0])
                error("%s", ctx.message);
        *N = ctx.N;
        *M = ctx.M;
}

void ancestrymap2geno(char *input_file, char *output_file, int *N, int *M)
{
        am_convert(input_file, output_file, N, M, 0);
}

void ancestrymap2lfmm(char *input_file, char *output_file, int *N, int *M)
{
        am_convert(input_file, output_file, N, M, 1);
}

int nb_ind_ancestrymap(char *input_file)
{
        am_context ctx = {0};
        ctx.input_name = input_file ? input_file : "<null>";
        ctx.input = input_file && *input_file ? fopen(input_file, "rb") : NULL;
        if (!ctx.input)
                am_fail(&ctx, "Unable to open ANCESTRYMAP input.");
        else
                am_count(&ctx, 1);
        am_cleanup(&ctx);
        if (ctx.message[0])
                error("%s", ctx.message);
        return ctx.N;
}

void read_ancestrymap(char *input_file, int N, int M, int *data)
{
        am_context ctx = {0};
        ctx.input_name = input_file ? input_file : "<null>";
        ctx.N = N;
        ctx.M = M;
        ctx.data = data; /* borrowed; cleanup must not free caller storage */
        if (!data || !am_matrix_size(&ctx)) {
                if (!ctx.message[0]) am_fail(&ctx, "Invalid ANCESTRYMAP data buffer.");
        } else {
                ctx.input = input_file && *input_file ? fopen(input_file, "rb") : NULL;
                if (!ctx.input) am_fail(&ctx, "Unable to open ANCESTRYMAP input.");
                else am_read_matrix(&ctx);
        }
        am_cleanup(&ctx);
        am_warning(&ctx);
        if (ctx.message[0]) error("%s", ctx.message);
}

void read_line_ancestrymap(char *szbuff, int *allele, char *name, int i, int j,
                           char *input, int *warning)
{
        am_context ctx = {0};
        char *token;
        int value;
        ctx.input_name = input ? input : "<null>";
        if (!allele || !name || !warning)
                am_fail(&ctx, "Invalid ANCESTRYMAP line arguments.");
        else if (am_parse(&ctx, szbuff, &token, &value)) {
                memcpy(name, token, strlen(token) + 1);
                *allele = value;
                if (!*warning && am_noncanonical(value)) {
                        ctx.warned = 1;
                        ctx.warning_i = (size_t)i;
                        ctx.warning_j = (size_t)j;
                        *warning = 1;
                }
        }
        if (ctx.message[0]) error("%s", ctx.message);
        am_warning(&ctx);
}
