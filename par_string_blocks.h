// STRING_BLOCKS :: https://github.com/prideout/par
// String extraction and concatenation, especially useful for snippets of GLSL or Lua.
//
// This little library extracts blocks of text from a memory blob or file, then lets you retrieve
// them by name or dump them out to a C header. It also makes it easy to glue together a sequence of
// blocks.
//
// Each block of text is assigned a name using a prefix line that starts with three dash characters,
// such as "--- the_name" or "--- my.block".
//
// For example, suppose you have a file called "shaders.glsl" that looks like this:
//
//    --- my_shader
//    void main() { ... }
//    --- common
//    uniform vec4 resolution;
//    uniform vec4 color;
//
// You can use this library to read in the file and extract one of the blocks:
//
//     parsb_context* blocks = parsb_create_context((parsb_options){});
//     parsb_add_blocks_from_file(blocks, "shaders.glsl");
//     const char* single = parsb_get_blocks(blocks, "my_shader");
//
// You can also concatenate blocks using a space-delimited list of names:
//
//     const char* concatenated = parsb_get_blocks(blocks, "common my_shader");
//
// You can also add or replace blocks on the fly:
//
//     parsb_add_block(blocks, "prefix", "#version 330\n");
//     const char* concatenated = parsb_get_blocks(blocks, "prefix common my_shader");
//
// The "blocks" context in the above examples holds a cache of generated strings, so be sure to
// destroy it when you're done:
//
//     parsb_destroy_context(blocks);
//
// Distributed under the MIT License, see bottom of file.

#ifndef PAR_STRING_BLOCKS_H
#define PAR_STRING_BLOCKS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

// OPTIONS
// -------
// line_directives ... adds #line annotations into concatenated strings for better error messages.
typedef struct parsb_options {
    bool line_directives;
} parsb_options;

// CONTEXT CREATION AND DESTRUCTION
// --------------------------------
// A context is an opaque handle to a memory arena. All generated strings are owned by the context
// and freed when the context is destroyed.
typedef struct parsb_context_s parsb_context;
parsb_context* parsb_create_context(parsb_options);
void parsb_destroy_context(parsb_context*);

// ADDING AND REPLACING BLOCKS
// ---------------------------
// When using the plural form (add_blocks), the submitted buffer may contain multiple blocks, each
// with a name defined by its closest preceding triple-dash line. If a block with the specified name
// already exists, it gets replaced.
//
// The singular form (add_block) adds a single block whose name is explicitly specified as an
// argument. Again, if a block with the given name already exists, it gets replaced.
//
// These functions do not retain the passed-in strings so clients can free them after pushing them.
void parsb_add_blocks(parsb_context*, const char* buffer, int buffer_size);
void parsb_add_block(parsb_context*, const char* name, const char* body);
#ifndef PARSB_NO_STDIO
void parsb_add_blocks_from_file(parsb_context* context, const char* filename);
#endif

// EXTRACTING AND CONCATENATING BLOCKS
// -----------------------------------
// The block_names string is a space-separated list of block names that are being requested. The
// returned string is owned by the context, so please make a copy if you need it to outlive the
// context. If the returned string is null, then one or more of the block names could not be found.
const char* parsb_get_blocks(parsb_context*, const char* block_names);

// GETTING BLOCKS BY INDEX
// -----------------------
int parsb_get_num_blocks(const parsb_context*);
void parsb_get_block(const parsb_context*, int index, const char** name, const char** body);

// SAVING THE BLOCK LIST
// ---------------------
// These functions export the entire "database" of atomic blocks.
typedef void (*parsb_write_line)(const char* line, void* userdata);
void parsb_write_blocks(parsb_context*, parsb_write_line writefn, void* user);

#ifndef PARSB_MAX_NUM_BLOCKS
#define PARSB_MAX_NUM_BLOCKS 128
#endif

#ifndef PARSB_MAX_NAME_LENGTH
#define PARSB_MAX_NAME_LENGTH 256
#endif

#ifndef PARSB_MAX_LINE_LENGTH
#define PARSB_MAX_LINE_LENGTH 256
#endif

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_STRING_BLOCKS_IMPLEMENTATION

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#ifndef PARSB_NO_STDIO
#include <stdio.h>
#endif

#define PARSB_MIN(a, b) (a > b ? b : a) 

typedef struct {
    int count;
    char* values[PARSB_MAX_NUM_BLOCKS];
    char* names[PARSB_MAX_NUM_BLOCKS];
} parsb__list;

struct parsb_context_s {
    parsb_options options;
    parsb__list blocks;
    parsb__list results;
};

static char* parsb__add_or_replace(parsb_context*, const char* id, const char* value,
    int value_size, int line_number);
static char* parsb__list_add(parsb__list*, const char* id, const char* value, int value_size,
    int line_number);
static char* parsb__list_get(parsb__list*, const char* id, int idlen);
static void parsb__list_free(parsb__list* );

parsb_context* parsb_create_context(parsb_options options) {
    parsb_context* context = (parsb_context*) calloc(1, sizeof(parsb_context));
    context->options = options;
    return context;
}

void parsb_destroy_context(parsb_context* context) {
    parsb__list_free(&context->blocks);
    parsb__list_free(&context->results);
    free(context);
}

void parsb_add_blocks(parsb_context* context, const char* blob, int buffer_size) {
    const char* previous_block = 0;
    char previous_name[PARSB_MAX_NAME_LENGTH];
    int line_number = 0;
    int block_line_number = 0;
    for (int i = 0; i < buffer_size - 4; i++) {
        if (blob[i] != '-' || blob[i + 1] != '-' || blob[i + 2] != '-' || blob[i + 3] != ' ') {
            if (blob[i] == '\n') {
                line_number++;
            }
            continue;
        }
        if (previous_block) {
            parsb__add_or_replace(context, previous_name, previous_block,
                i - (previous_block - blob), block_line_number);
        }
        i += 4;
        const char* name = blob + i;
        const char* block_start = 0;
        for (; i < buffer_size; i++) {
            if (blob[i] == '\n') {
                line_number++;
                int name_length = i - (name - blob);
                memcpy(previous_name, name, name_length);
                block_line_number = line_number + 2;
                previous_name[name_length] = 0;
                block_start = blob + i + 1;
                break;
            }
            if (isspace(blob[i])) {
                int name_length = i - (name - blob);
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
        parsb__add_or_replace(context, previous_name, previous_block,
            buffer_size - (previous_block - blob), block_line_number);
    }
}

void parsb_add_block(parsb_context* context, const char* name, const char* body) {
    char* dup = strdup(body);
    parsb__add_or_replace(context, name, dup, 1 + strlen(body), 0);
}

const char* parsb_get_blocks(parsb_context* context, const char* block_names) {
    int len = strlen(block_names);
    const char* name = block_names;
    int name_length = 0;
    int result_length = 0;

    // First pass determines the amount of required memory.
    int num_names = 0;
    for (int i = 0; i < len; i++) {
        char c = block_names[i];
        if (isspace(c) || !c) {
            const char* block = parsb__list_get(&context->blocks, name, name_length);
            if (block) {
                result_length += strlen(block);
                num_names++;
            } else {
                return NULL;
            }
            name_length = 0;
            name = block_names + i + 1;
        } else {
            name_length++;
        }
    }
    const char* block = parsb__list_get(&context->blocks, name, name_length);
    if (block) {
        result_length += strlen(block);
        num_names++;
    }

    // If no concatenation is required, return early.
    if (num_names == 1) {
        return parsb__list_get(&context->blocks, name, name_length);
    }

    // Allocate storage for the result.
    char* result = parsb__list_add(&context->results, 0, 0, result_length, 0);
    char* cursor = result;

    // Second pass populates the result.
    name = block_names;
    name_length = 0;
    for (int i = 0; i < len; i++) {
        char c = block_names[i];
        if (isspace(c) || !c) {
            const char* block = parsb__list_get(&context->blocks, name, name_length);
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
    block = parsb__list_get(&context->blocks, name, name_length);
    if (block) {
        memcpy(cursor, block, strlen(block));
        cursor += strlen(block);
    }
    return result;
}

int parsb_get_num_blocks(const parsb_context* context) {
    return context->blocks.count;
}

void parsb_get_block(const parsb_context* context, int index, const char** name,
    const char** body) {
    if (index < 0 || index >= context->blocks.count) {
        return;
    }
    *name = context->blocks.names[index];
    *body = context->blocks.values[index];
}

void parsb_write_blocks(parsb_context* context, parsb_write_line writefn, void* userdata) {
    char line[PARSB_MAX_LINE_LENGTH + 1] = {0};
    for (int i = 0; i < context->blocks.count; i++) {

        sprintf(line, "--- %s", context->blocks.names[i]);
        writefn(line, userdata);

        const char* cursor = context->blocks.values[i];
        const int blocklen = strlen(cursor);
        int previous = 0;
        for (int i = 0; i < blocklen; i++) {
            if (cursor[i] == '\n') {
                int line_length = PARSB_MIN(i - previous, PARSB_MAX_LINE_LENGTH);
                memcpy(line, cursor + previous, line_length);
                line[line_length] = 0;
                writefn(line, userdata);
                previous = i + 1;
            } else if (i == blocklen - 1) {
                int line_length = PARSB_MIN(1 + i - previous, PARSB_MAX_LINE_LENGTH);
                memcpy(line, cursor + previous, line_length);
                line[line_length] = 0;
                writefn(line, userdata);
                previous = i + 1;
            }
        }
    }
}

static char* parsb__add_or_replace(parsb_context* context, const char* id, const char* value,
    int value_size, int line_number) {
    line_number = context->options.line_directives ? line_number : 0;
    const size_t idlen = strlen(id);
    for (int i = 0; i < context->blocks.count; i++) {
        if (strncmp(id, context->blocks.names[i], idlen) == 0) {
            free(context->blocks.values[i]);
            context->blocks.values[i] = strndup(value, value_size);
            return context->blocks.values[i];
        }
    }
    return parsb__list_add(&context->blocks, id, value, value_size, line_number);
}

static char* parsb__list_add(parsb__list* list, const char* name,
    const char* value, int value_size, int line_number) {
    if (value_size == 0) {
        return NULL;
    }

    if (list->count == PARSB_MAX_NUM_BLOCKS) {
        assert(false && "Please increase PARSB_MAX_NUM_BLOCKS.");
        return NULL;
    }

    char* storage;
    char* cursor;

    if (line_number > 0) {
        char line_directive[16] = {0};
        int prefix_length = snprintf(line_directive, 16, "\n#line %d\n", line_number);
        storage = (char*) calloc(1, prefix_length + value_size + 1);
        memcpy(storage, line_directive, prefix_length);
        cursor = storage + prefix_length;
    } else {
        storage = cursor = (char*) calloc(1, value_size + 1);
    }

    if (value) {
        memcpy(cursor, value, value_size);
    }

    #if PARSB_ENABLE_TRIM
    value_size--;
    while (isspace(cursor[value_size])) {
        cursor[value_size] = 0;
        value_size--;
        if (value_size == 0) {
            break;
        }
    }
    #endif

    if (name) {
        list->names[list->count] = strdup(name);
    } else {
        list->names[list->count] = 0;
    }

    list->values[list->count] = storage;
    list->count++;
    return storage;
}

static char* parsb__list_get(parsb__list* list, const char* name, int idlen) {
    for (int i = 0; i < list->count; i++) {
        if (strncmp(name, list->names[i], idlen) == 0) {
            return list->values[i];
        }
    }
    return NULL;
}

static void parsb__list_free(parsb__list* list) {
    for (int i = 0; i < list->count; i++) {
        free(list->names[i]);
        free(list->values[i]);
    }
    list->count = 0;
}

#ifndef PARSB_NO_STDIO

void parsb_add_blocks_from_file(parsb_context* context, const char* filename) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Unable to open %s\n", filename);
        return;
    }
    fseek(f, 0, SEEK_END);
    int length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buffer = (char*) malloc(length);
    fread(buffer, 1, length, f);
    fclose(f);
    parsb_add_blocks(context, buffer, length);
    free(buffer);
}

#endif

#endif // PAR_STRING_BLOCKS_IMPLEMENTATION
#endif // PAR_STRING_BLOCKS_H

// par_string_blocks is distributed under the MIT license:
//
// Copyright (c) 2020 Philip Rideout
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
