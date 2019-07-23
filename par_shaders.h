// SHADERS :: https://github.com/prideout/par
// String extraction and concatenation utilities designed for shaders.
//
// This little library extracts blocks of text from a memory blob,
// then lets you retrieve them by name or dump them out to a C header.
// It also makes it easy to glue together a sequence of blocks.
//
// Each block of text is assigned a name using a prefix line that starts with
// two dash characters, such as "-- the_name" or "-- my.block".
//
// For example, the contents of your memory blob could look like this:
//
//   -- my_vs
//   void main() { ... }
//
//   -- my_decls
//   uniform vec4 resolution;
//   uniform vec4 color;
//
// You could then extract text blocks from the memory blob like so:
//
//   parsh_context* ctx = parsh_create_context({});
//   parsh_add_blocks(ctx, pointer_to_blob, size_of_blob_in_bytes);
//   free(pointer_to_blob);
//   parsh_add_block(ctx, "prefix", "#version 330\n");
//   const char* concatenated = parsh_get_blocks(ctx, "prefix my_decls my_vs");
//   ...use "concatenated" here...
//   parsh_destroy_context(ctx);
//
// This library is similar to the string wrangling library described
// in the following old post, but has been updated to be a no-dependency
// single-file library in the style of the STB libraries.
//
//   https://prideout.net/blog/old/blog/index.html@p=11.html
//
// Distributed under the MIT License, see bottom of file.

#ifndef PAR_SHADERS_H
#define PAR_SHADERS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
    bool enable_line_directives;
} parsh_config;

// Opaque handle to a memory arena. All generated strings are owned by the
// library, and freed when the context is destroyed.
typedef struct parsh_context_s parsh_context;

parsh_context* parsh_create_context(parsh_config);
void parsh_destroy_context(parsh_context*);
void parsh_add_blocks(parsh_context*, const char* buffer, size_t buffer_size);
void parsh_add_block(parsh_context*, const char* name, const char* body);
const char* parsh_get_blocks(parsh_context*, const char* block_names);

typedef void (*parsh_write_line)(const char* line, void* userdata);
void parsh_write_cstring(parsh_context*, parsh_write_line writefn, void* user);

parsh_context* parsh_create_context_from_file(const char* filename);
void parsh_add_blocks_from_file(parsh_context* context, const char* filename);

#ifndef PARSH_MAX_NUM_BLOCKS
#define PARSH_MAX_NUM_BLOCKS 128
#endif

#ifndef PARSH_MAX_NAME_LENGTH
#define PARSH_MAX_NAME_LENGTH 256
#endif

#ifndef PARSH_MAX_LINE_LENGTH
#define PARSH_MAX_LINE_LENGTH 256
#endif

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_SHADERS_IMPLEMENTATION

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PAR_PI
#define PAR_PI (3.14159265359)
#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#define PAR_SQR(a) ((a) * (a))
#endif

typedef struct {
    size_t count;
    char* values[PARSH_MAX_NUM_BLOCKS];
    char* names[PARSH_MAX_NUM_BLOCKS];
} parsh__list;

struct parsh_context_s {
    parsh_config config;
    parsh__list blocks;
    parsh__list results;
};

static char* parsh__list_add(parsh__list*, const char* id,
    const char* value, size_t value_size, int line_number);
static const char* parsh__list_get(parsh__list*, const char* id, size_t idlen);
static void parsh__list_free(parsh__list* );

parsh_context* parsh_create_context(parsh_config config) {
    parsh_context* context = (parsh_context*) calloc(1, sizeof(parsh_context));
    context->config = config;
    return context;
}

void parsh_destroy_context(parsh_context* context) {
    parsh__list_free(&context->blocks);
    parsh__list_free(&context->results);
    free(context);
}

void parsh_add_blocks(parsh_context* context, const char* blob,
    size_t buffer_size) {
    const char* previous_block = 0;
    char previous_name[PARSH_MAX_NAME_LENGTH];
    int line_number = 0;
    int block_line_number = 0;
    for (size_t i = 0; i < buffer_size - 3; i++) {
        if (blob[i] != '-' || blob[i + 1] != '-' || blob[i + 2] != ' ') {
            if (blob[i] == '\n') {
                line_number++;
            }
            continue;
        }
        if (previous_block) {
            parsh__list_add(&context->blocks, previous_name, previous_block,
                i - (previous_block - blob),
                context->config.enable_line_directives ? block_line_number : 0);
        }
        i += 3;
        const char* name = blob + i;
        const char* block_start = 0;
        for (; i < buffer_size; i++) {
            if (blob[i] == '\n') {
                line_number++;
                size_t name_length = i - (name - blob);
                memcpy(previous_name, name, name_length);
                block_line_number = line_number + 2;
                previous_name[name_length] = 0;
                block_start = blob + i + 1;
                break;
            }
            if (isspace(blob[i])) {
                size_t name_length = i - (name - blob);
                memcpy(previous_name, name, name_length);
                block_line_number = line_number + 2;
                previous_name[name_length] = 0;
                for (i++; i < buffer_size; i++) {
                    if (blob[i] == '\n') {
                        line_number++;
                        block_start = blob + i + 1;
                        break;
                    }
                }
                break;
            }
        }
        if (block_start == 0) {
            return;
        }
        previous_block = block_start;
    }
    if (previous_block) {
        parsh__list_add(&context->blocks, previous_name, previous_block,
            buffer_size - (previous_block - blob),
            context->config.enable_line_directives ? block_line_number : 0);
    }
}

void parsh_add_block(parsh_context* context, const char* name,
    const char* body) {
    char* dup = (char*) malloc(strlen(body) + 1);
    memcpy(dup, body, 1 + strlen(body));
    parsh__list_add(&context->blocks, name, dup, 1 + strlen(body), 0);
}

const char* parsh_get_blocks(parsh_context* context, const char* block_names) {
    size_t len = strlen(block_names);
    const char* name = block_names;
    size_t name_length = 0;
    size_t result_length = 0;

    // First pass determines the amount of required memory.
    size_t num_names = 0;
    for (size_t i = 0; i < len; i++) {
        char c = block_names[i];
        if (isspace(c) || !c) {
            const char* block = parsh__list_get(&context->blocks, name,
                name_length);
            if (block) {
                result_length += strlen(block);
                num_names++;
            }
            name_length = 0;
            name = block_names + i + 1;
        } else {
            name_length++;
        }
    }
    const char* block = parsh__list_get(&context->blocks, name,
        name_length);
    if (block) {
        result_length += strlen(block);
        num_names++;
    }

    // If no concatenation is required, return early.
    if (num_names == 1) {
        return parsh__list_get(&context->blocks, name, name_length);
    }

    // Allocate storage for the result.
    char* result = parsh__list_add(&context->results, 0, 0, result_length, 0);
    char* cursor = result;

    // Second pass populates the result.
    name = block_names;
    name_length = 0;
    for (size_t i = 0; i < len; i++) {
        char c = block_names[i];
        if (isspace(c) || !c) {
            const char* block = parsh__list_get(&context->blocks, name,
                name_length);
            if (block) {
                memcpy(cursor, block, strlen(block));
                cursor += strlen(block);
            }
            name_length = 0;
            name = block_names + i + 1;
        } else {
            name_length++;
        }
    }
    block = parsh__list_get(&context->blocks, name, name_length);
    if (block) {
        memcpy(cursor, block, strlen(block));
        cursor += strlen(block);
    }
    return result;
}

void parsh_write_cstring(parsh_context* context, parsh_write_line writefn,
    void* userdata) {
    char line[PARSH_MAX_LINE_LENGTH + 4] = {0};
    for (size_t i = 0; i < context->blocks.count; i++) {

        sprintf(line, "\"-- %s\\n\"", context->blocks.names[i]);
        writefn(line, userdata);

        const char* cursor = context->blocks.values[i];
        const size_t blocklen = strlen(cursor);
        size_t previous = 0;
        for (size_t i = 0; i < blocklen; i++) {
            if (cursor[i] == '\n' || i == blocklen - 1) {
                size_t line_length = PAR_MIN(i - previous,
                    PARSH_MAX_LINE_LENGTH);
                if (i == blocklen - 1) {
                    line_length++;
                }
                line[0] = '\"';
                memcpy(line + 1, cursor + previous, line_length);
                line[1 + line_length] = '\\';
                line[2 + line_length] = 'n';
                line[3 + line_length] = '\"';
                line[4 + line_length] = 0;
                writefn(line, userdata);
                previous = i + 1;
            }
        }
    }
}

static char* parsh__list_add(parsh__list* list, const char* name,
    const char* value, size_t value_size, int line_number) {
    if (value_size == 0) {
        return 0;
    }
    if (list->count == PARSH_MAX_NUM_BLOCKS) {
        assert(false && "Please increase PARSH_MAX_NUM_BLOCKS.");
        return 0;
    }

    char* storage;
    char* cursor;

    if (line_number > 0) {
        char line_directive[16] = {0};
        size_t prefix_length =
            snprintf(line_directive, 16, "\n#line %d\n", line_number);
        storage = (char*) calloc(1, prefix_length + value_size + 1);
        memcpy(storage, line_directive, prefix_length);
        cursor = storage + prefix_length;
    } else {
        storage = cursor = (char*) calloc(1, value_size + 1);
    }

    if (value) {
        memcpy(cursor, value, value_size--);
    }

    while (isspace(cursor[value_size])) {
        cursor[value_size] = 0;
        value_size--;
        if (value_size == 0) {
            break;
        }
    }

    if (name) {
        char* dup = (char*) malloc(strlen(name) + 1);
        memcpy(dup, name, strlen(name) + 1);
        list->names[list->count] = dup;
    } else {
        list->names[list->count] = 0;
    }

    list->values[list->count] = storage;
    list->count++;
    return storage;
}

static const char* parsh__list_get(parsh__list* list, const char* name,
    size_t idlen) {
    for (size_t i = 0; i < list->count; i++) {
        if (strncmp(name, list->names[i], idlen) == 0) {
            return list->values[i];
        }
    }
    return 0;
}

static void parsh__list_free(parsh__list* list) {
    for (size_t i = 0; i < list->count; i++) {
        free(list->names[i]);
        free(list->values[i]);
    }
    list->count = 0;
}

#ifdef PARSH_ENABLE_STDIO

parsh_context* parsh_create_context_from_file(const char* filename) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        return NULL;
    }
    fseek(f, 0, SEEK_END);
    size_t length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buffer = malloc(length);
    fread(buffer, 1, length, f);
    fclose(f);
    parsh_context* shaders = parsh_create_context((parsh_config){
        .enable_line_directives = true
    });
    parsh_add_blocks(shaders, buffer, length);
    free(buffer);
    return shaders;
}

void parsh_add_blocks_from_file(parsh_context* context, const char* filename) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Unable to open %s\n", filename);
        return;
    }
    fseek(f, 0, SEEK_END);
    size_t length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buffer = malloc(length);
    fread(buffer, 1, length, f);
    fclose(f);
    parsh_add_blocks(context, buffer, length);
    free(buffer);
}

#endif

#ifdef PARSH_ENABLE_MAIN

void write_line(const char* ln, void* userdata) {
    FILE* outfile = (FILE*) userdata;
    fputs(ln, outfile);
    fputc('\n', outfile);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        puts("Usage: parsh srcfile dstfile array_name");
        return 1;
    }
    const char* srcfile = argv[1];
    const char* dstfile = argv[2];
    const char* array_name = argv[3];

    FILE *f = fopen(srcfile, "rb");
    fseek(f, 0, SEEK_END);
    size_t length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buffer = malloc (length);
    fread(buffer, 1, length, f);
    fclose(f);

    parsh_context* ctx = parsh_create_context((parsh_config){
        .enable_line_directives = true
    });
    parsh_add_blocks(ctx, buffer, length);
    free(buffer);

    FILE* outfile = fopen(dstfile, "wt");
    fprintf(outfile, "const char %s[] = \n", array_name);
    parsh_write_cstring(ctx, write_line, outfile);
    fprintf(outfile, ";\n");
    fclose(outfile);
    parsh_destroy_context(ctx);
    return 0;
}

#endif

#endif // PAR_SHADERS_IMPLEMENTATION
#endif // PAR_SHADERS_H

// par_shaders is distributed under the MIT license:
//
// Copyright (c) 2019 Philip Rideout
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
