// STREAMLINES :: https://github.com/prideout/par
// Simple C library for triangulating wide lines, curves, and streamlines.
//
// Documentation at https://prideout.net/blog/par_streamlines/
//
// The MIT License
// Copyright (c) 2019 Philip Rideout

#ifndef PAR_STREAMLINES_H
#define PAR_STREAMLINES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

typedef struct {
    float x;
    float y;
} par_streamlines_position;

typedef struct {
    float distance_along_spine;       // non-normalized distance along the entire curve
    float line_length;                // the length of the entire curve
    float signed_distance_from_spine; // either + or - depending on the side
    uint16_t line_index;              // tells which curve this vertex belongs to
    uint16_t segment_index;           // vertices are shared so this might be interpolated
} par_streamlines_annotation;

typedef struct {
    uint32_t num_vertices;
    uint32_t num_triangles;
    par_streamlines_position* vertex_positions;
    par_streamlines_annotation* vertex_annotations;
    uint32_t* triangle_indices;
} par_streamlines_mesh;

typedef struct {
    float lines_thickness;
    uint32_t curves_level_of_detail;
    float streamlines_seed_spacing;
    float streamlines_seed_viewport[4];
    uint32_t streamlines_num_frames;
} par_streamlines_config;

typedef struct {
    uint32_t num_vertices;
    uint32_t num_lines;
    par_streamlines_position* vertices;
    uint32_t* line_lengths;
} par_streamlines_path_list;

typedef struct par_streamlines_context_s par_streamlines_context;

typedef void (*par_streamlines_callback)(float domain[2], float range[2]);

par_streamlines_context* par_streamlines_create_context(par_streamlines_config config);
void par_streamlines_destroy_context(par_streamlines_context* context);
par_streamlines_mesh* par_streamlines_draw_lines(par_streamlines_context* context,
    par_streamlines_path_list paths);
par_streamlines_mesh* par_streamlines_draw_curves_cubic(par_streamlines_context* context,
    par_streamlines_path_list paths);
par_streamlines_mesh* par_streamlines_draw_curves_quadratic(par_streamlines_context* context,
    par_streamlines_path_list paths);
par_streamlines_mesh* par_streamlines_draw_streamlines(par_streamlines_context* context,
    par_streamlines_callback func, uint32_t frame_index);

#ifdef __cplusplus
}
#endif

#ifndef PAR_MALLOC
#define PAR_MALLOC(T, N) ((T*) malloc(N * sizeof(T)))
#define PAR_CALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))
#define PAR_REALLOC(T, BUF, N) ((T*) realloc(BUF, sizeof(T) * (N)))
#define PAR_FREE(BUF) free(BUF)
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_STREAMLINES_IMPLEMENTATION

#include <stdlib.h>

struct par_streamlines_context_s {
    par_streamlines_config config;
    par_streamlines_mesh result;
};

par_streamlines_context* par_streamlines_create_context(par_streamlines_config config)
{
    par_streamlines_context* context = PAR_CALLOC(par_streamlines_context, 1);
    context->config = config;
    printf("LEFT = %f\n", context->config.streamlines_seed_viewport[0]);
    return context;
}

void par_streamlines_destroy_context(par_streamlines_context* context)
{
    PAR_FREE(context);
}

par_streamlines_mesh* par_streamlines_draw_lines(par_streamlines_context* context,
    par_streamlines_path_list paths)
{
    return &context->result;
}

par_streamlines_mesh* par_streamlines_draw_curves_cubic(par_streamlines_context* context,
    par_streamlines_path_list paths)
{
    return &context->result;
}

par_streamlines_mesh* par_streamlines_draw_curves_quadratic(par_streamlines_context* context,
    par_streamlines_path_list paths)
{
    return &context->result;
}

par_streamlines_mesh* par_streamlines_draw_streamlines(par_streamlines_context* context,
    par_streamlines_callback func, uint32_t frame_index)
{
    return &context->result;
}

#endif // PAR_STREAMLINES_IMPLEMENTATION
#endif // PAR_STREAMLINES_H
