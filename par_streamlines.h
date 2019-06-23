// STREAMLINES :: https://github.com/prideout/par
// Simple C library for triangulating wide lines, curves, and streamlines.
//
// Documentation at https://prideout.net/blog/par_streamlines/
// Simple example:
//
//   #define PAR_STREAMLINES_IMPLEMENTATION
//   #include "par_streamlines.h"
//
//   par_streamlines_config cfg = { .thickness = 3 };
//   par_streamlines_context* ctx = par_streamlines_create_context(cfg);
//   par_streamlines_position vertices[] = { {0, 0}, {2, 1}, {4, 0} };
//   uint16_t spine_lengths[] = { 3 };
//   par_streamlines_spine_list spines = {
//       .num_vertices = 3,
//       .num_spines = 1,
//       .vertices = vertices,
//       .spine_lengths = spine_lengths
//   };
//   par_streamlines_mesh* mesh = par_streamlines_draw_lines(ctx, spines);
//   ...
//   par_streamlines_destroy_context(ctx);
//
// The MIT License
// Copyright (c) 2019 Philip Rideout

#ifndef PAR_STREAMLINES_H
#define PAR_STREAMLINES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>

typedef struct {
    float x;
    float y;
} par_streamlines_position;

// TODO: change this to u v direction_x direction_y
typedef struct {
    float distance_along_spine;       // non-normalized distance along the curve
    float spine_length;               // the length of the entire curve
    float signed_distance_from_spine; // either + or - depending on the side
    uint16_t spine_index;             // tells which curve this vertex lives in
    uint16_t segment_index;           // note that this might be interpolated
} par_streamlines_annotation;

typedef struct {
    uint32_t num_vertices;
    uint32_t num_triangles;
    par_streamlines_position* vertex_positions;
    par_streamlines_annotation* vertex_annotations;
    uint32_t* triangle_indices;
} par_streamlines_mesh;

typedef struct {
    float thickness;
    uint32_t curves_level_of_detail;
    float streamlines_seed_spacing;
    float streamlines_seed_viewport[4];
    uint32_t streamlines_num_frames;
    bool wireframe; // creates 4 indices per triangle instead of 3
} par_streamlines_config;

typedef struct {
    uint32_t num_vertices;
    uint16_t num_spines;
    par_streamlines_position* vertices;
    uint16_t* spine_lengths;
    bool closed;
} par_streamlines_spine_list;

typedef struct par_streamlines_context_s par_streamlines_context;

typedef void (*par_streamlines_callback)(float domain[2], float range[2]);

par_streamlines_context* par_streamlines_create_context(
    par_streamlines_config config);

void par_streamlines_destroy_context(par_streamlines_context* context);

par_streamlines_mesh* par_streamlines_draw_lines(
    par_streamlines_context* context, par_streamlines_spine_list spines);

par_streamlines_mesh* par_streamlines_draw_curves_cubic(
    par_streamlines_context* context, par_streamlines_spine_list spines);

par_streamlines_mesh* par_streamlines_draw_curves_quadratic(
    par_streamlines_context* context, par_streamlines_spine_list spines);

par_streamlines_mesh* par_streamlines_draw_streamlines(
    par_streamlines_context* context, par_streamlines_callback func,
    uint32_t frame_index);

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_STREAMLINES_IMPLEMENTATION

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifndef PAR_PI
#define PAR_PI (3.14159265359)
#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#define PAR_SQR(a) ((a) * (a))
#endif

#ifndef PAR_MALLOC
#define PAR_MALLOC(T, N) ((T*) malloc(N * sizeof(T)))
#define PAR_CALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))
#define PAR_REALLOC(T, BUF, N) ((T*) realloc(BUF, sizeof(T) * (N)))
#define PAR_FREE(BUF) free(BUF)
#endif

#ifndef PAR_ARRAY
#define pa_free(a) ((a) ? PAR_FREE(pa___raw(a)), 0 : 0)
#define pa_push(a, v) (pa___maybegrow(a, 1), (a)[pa___n(a)++] = (v))
#define pa_count(a) ((a) ? pa___n(a) : 0)
#define pa_add(a, n) (pa___maybegrow(a, n), pa___n(a) += (n))
#define pa_last(a) ((a)[pa___n(a) - 1])
#define pa_end(a) (a + pa_count(a))
#define pa_clear(arr) if (arr) pa___n(arr) = 0

#define pa___raw(a) ((int*) (a) -2)
#define pa___m(a) pa___raw(a)[0]
#define pa___n(a) pa___raw(a)[1]
#define pa___needgrow(a, n) ((a) == 0 || pa___n(a) + (n) >= pa___m(a))
#define pa___maybegrow(a, n) (pa___needgrow(a, (n)) ? pa___grow(a, n) : 0)
#define pa___grow(a, n) (*((void**)& (a)) = pa___growf((void*) (a), (n), \
        sizeof(*(a))))

static void* pa___growf(void* arr, int increment, int itemsize)
{
    int dbl_cur = arr ? 2 * pa___m(arr) : 0;
    int min_needed = pa_count(arr) + increment;
    int m = dbl_cur > min_needed ? dbl_cur : min_needed;
    int* p = PAR_REALLOC(int, arr ? pa___raw(arr) : 0,
        itemsize * m / sizeof(int) + 2);
    if (p) {
        if (!arr) {
            p[1] = 0;
        }
        p[0] = m;
        return p + 2;
    }
    return (void*) (2 * sizeof(int));
}

#endif

struct par_streamlines_context_s {
    par_streamlines_config config;
    par_streamlines_mesh result;
};

par_streamlines_context* par_streamlines_create_context(
    par_streamlines_config config)
{
    par_streamlines_context* context = PAR_CALLOC(par_streamlines_context, 1);
    context->config = config;
    return context;
}

void par_streamlines_destroy_context(par_streamlines_context* context)
{
    pa_free(context->result.vertex_annotations);
    pa_free(context->result.vertex_positions);
    PAR_FREE(context);
}

par_streamlines_mesh* par_streamlines_draw_lines(
    par_streamlines_context* context, par_streamlines_spine_list spines)
{
    typedef par_streamlines_position Position;
    typedef par_streamlines_annotation Annotation;

    par_streamlines_mesh* mesh = &context->result;
    const float extent = context->config.thickness / 2.0f;
    const bool closed = spines.closed;
    const bool wireframe = context->config.wireframe;
    const uint32_t ind_per_tri = wireframe ? 4 : 3;

    mesh->num_vertices = 0;
    mesh->num_triangles = 0;

    for (uint32_t spine = 0; spine < spines.num_spines; spine++) {
        assert(spines.spine_lengths[spine] > 1);
        mesh->num_vertices += 2 * spines.spine_lengths[spine];
        mesh->num_triangles += 2 * (spines.spine_lengths[spine] - 1);
        if (closed) {
            mesh->num_vertices += 2;
            mesh->num_triangles += 2;
        }
    }

    pa_clear(mesh->vertex_annotations);
    pa_clear(mesh->vertex_positions);
    pa_clear(mesh->triangle_indices);

    pa_add(mesh->vertex_annotations, mesh->num_vertices);
    pa_add(mesh->vertex_positions, mesh->num_vertices);
    pa_add(mesh->triangle_indices, ind_per_tri * mesh->num_triangles);

    Annotation* dst_annotations = mesh->vertex_annotations;
    Position* dst_positions = mesh->vertex_positions;
    uint32_t* dst_indices = mesh->triangle_indices;

    const Position* src_position = spines.vertices;
    uint32_t base_index = 0;

    for (uint16_t spine = 0; spine < spines.num_spines; spine++) {
        const uint16_t spine_length = spines.spine_lengths[spine];
        const float dx = src_position[1].x - src_position[0].x;
        const float dy = src_position[1].y - src_position[0].y;
        const float segment_length = sqrtf(dx * dx + dy * dy);
        const float nx = -dy / segment_length;
        const float ny = dx / segment_length;

        float pnx = nx;
        float pny = ny;

        const Position first_src_position = src_position[0];
        const Position last_src_position = src_position[spine_length - 1];

        if (closed) {
            const float dx = src_position[0].x - last_src_position.x;
            const float dy = src_position[0].y - last_src_position.y;
            const float segment_length = sqrtf(dx * dx + dy * dy);
            const float nx = -dy / segment_length;
            const float ny = dx / segment_length;

            pnx += nx;
            pny += ny;
            const float invlen = 1.0f / sqrtf(pnx * pnx + pny * pny);
            pnx *= invlen;
            pny *= invlen;
        }

        dst_positions[0].x = src_position[0].x + pnx * extent;
        dst_positions[0].y = src_position[0].y + pny * extent;
        dst_positions[1].x = src_position[0].x - pnx * extent;
        dst_positions[1].y = src_position[0].y - pny * extent;

        const Position first_dst_positions[2] = {
            dst_positions[0],
            dst_positions[1]
        };

        src_position++;
        dst_positions += 2;

        dst_annotations[0].distance_along_spine = 0;
        dst_annotations[0].signed_distance_from_spine = extent;
        dst_annotations[0].spine_index = (uint16_t) spine;
        dst_annotations[0].segment_index = 0;
        dst_annotations[1] = dst_annotations[0];
        dst_annotations[1].signed_distance_from_spine = -extent;
        dst_annotations += 2;

        float distance_along_spine = segment_length;

        uint16_t segment_index = 1;
        for (; segment_index < spine_length - 1; segment_index++) {

            const float dx = src_position[1].x - src_position[0].x;
            const float dy = src_position[1].y - src_position[0].y;
            const float segment_length = sqrtf(dx * dx + dy * dy);
            const float nx = -dy / segment_length;
            const float ny = dx / segment_length;

            pnx += nx;
            pny += ny;
            const float invlen = 1.0f / sqrtf(pnx * pnx + pny * pny);
            pnx *= invlen;
            pny *= invlen;

            dst_positions[0].x = src_position[0].x + pnx * extent;
            dst_positions[0].y = src_position[0].y + pny * extent;
            dst_positions[1].x = src_position[0].x - pnx * extent;
            dst_positions[1].y = src_position[0].y - pny * extent;
            src_position++;
            dst_positions += 2;

            dst_annotations[0].distance_along_spine = distance_along_spine;
            dst_annotations[0].signed_distance_from_spine = extent;
            dst_annotations[0].spine_index = (uint16_t) spine;
            dst_annotations[0].segment_index = segment_index;
            dst_annotations[1] = dst_annotations[0];
            dst_annotations[1].signed_distance_from_spine = -extent;
            dst_annotations += 2;
            distance_along_spine += segment_length;

            if (wireframe) {
                dst_indices[0] = base_index + (segment_index - 1) * 2;
                dst_indices[1] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[2] = base_index + (segment_index - 0) * 2;
                dst_indices[3] = base_index + (segment_index - 1) * 2;

                dst_indices[4] = base_index + (segment_index - 0) * 2;
                dst_indices[5] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[6] = base_index + (segment_index - 0) * 2 + 1;
                dst_indices[7] = base_index + (segment_index - 0) * 2;
                dst_indices += 8;
            } else {
                dst_indices[0] = base_index + (segment_index - 1) * 2;
                dst_indices[1] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[2] = base_index + (segment_index - 0) * 2;

                dst_indices[3] = base_index + (segment_index - 0) * 2;
                dst_indices[4] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[5] = base_index + (segment_index - 0) * 2 + 1;
                dst_indices += 6;
            }
        }

        if (closed) {
            const float dx = first_src_position.x - src_position[0].x;
            const float dy = first_src_position.y - src_position[0].y;
            const float segment_length = sqrtf(dx * dx + dy * dy);
            const float nx = -dy / segment_length;
            const float ny = dx / segment_length;

            pnx += nx;
            pny += ny;
            const float invlen = 1.0f / sqrtf(pnx * pnx + pny * pny);
            pnx *= invlen;
            pny *= invlen;
        }

        dst_positions[0].x = src_position[0].x + pnx * extent;
        dst_positions[0].y = src_position[0].y + pny * extent;
        dst_positions[1].x = src_position[0].x - pnx * extent;
        dst_positions[1].y = src_position[0].y - pny * extent;
        src_position++;
        dst_positions += 2;

        dst_annotations[0].distance_along_spine = distance_along_spine;
        dst_annotations[0].signed_distance_from_spine = extent;
        dst_annotations[0].spine_index = (uint16_t) spine;
        dst_annotations[0].segment_index = segment_index;
        dst_annotations[1] = dst_annotations[0];
        dst_annotations[1].signed_distance_from_spine = -extent;
        dst_annotations += 2;
        distance_along_spine += segment_length;

        if (wireframe) {
            dst_indices[0] = base_index + (segment_index - 1) * 2;
            dst_indices[1] = base_index + (segment_index - 1) * 2 + 1;
            dst_indices[2] = base_index + (segment_index - 0) * 2;
            dst_indices[3] = base_index + (segment_index - 1) * 2;

            dst_indices[4] = base_index + (segment_index - 0) * 2;
            dst_indices[5] = base_index + (segment_index - 1) * 2 + 1;
            dst_indices[6] = base_index + (segment_index - 0) * 2 + 1;
            dst_indices[7] = base_index + (segment_index - 0) * 2;
            dst_indices += 8;
        } else {
            dst_indices[0] = base_index + (segment_index - 1) * 2;
            dst_indices[1] = base_index + (segment_index - 1) * 2 + 1;
            dst_indices[2] = base_index + (segment_index - 0) * 2;

            dst_indices[3] = base_index + (segment_index - 0) * 2;
            dst_indices[4] = base_index + (segment_index - 1) * 2 + 1;
            dst_indices[5] = base_index + (segment_index - 0) * 2 + 1;
            dst_indices += 6;
        }

        if (closed) {
            segment_index++;

            dst_positions[0] = first_dst_positions[0];
            dst_positions[1] = first_dst_positions[1];
            dst_positions += 2;

            dst_annotations[0].distance_along_spine = distance_along_spine;
            dst_annotations[0].signed_distance_from_spine = extent;
            dst_annotations[0].spine_index = (uint16_t) spine;
            dst_annotations[0].segment_index = segment_index;
            dst_annotations[1] = dst_annotations[0];
            dst_annotations[1].signed_distance_from_spine = -extent;
            dst_annotations += 2;
            distance_along_spine += segment_length;

            if (wireframe) {
                dst_indices[0] = base_index + (segment_index - 1) * 2;
                dst_indices[1] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[2] = base_index + (segment_index - 0) * 2;
                dst_indices[3] = base_index + (segment_index - 1) * 2;

                dst_indices[4] = base_index + (segment_index - 0) * 2;
                dst_indices[5] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[6] = base_index + (segment_index - 0) * 2 + 1;
                dst_indices[7] = base_index + (segment_index - 0) * 2;
                dst_indices += 8;
            } else {
                dst_indices[0] = base_index + (segment_index - 1) * 2;
                dst_indices[1] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[2] = base_index + (segment_index - 0) * 2;

                dst_indices[3] = base_index + (segment_index - 0) * 2;
                dst_indices[4] = base_index + (segment_index - 1) * 2 + 1;
                dst_indices[5] = base_index + (segment_index - 0) * 2 + 1;
                dst_indices += 6;
            }

            dst_annotations -= 2 + spine_length * 2;
            for (uint16_t i = 0; i < spine_length + 1; i++) {
                dst_annotations[0].spine_length = distance_along_spine;
                dst_annotations[1].spine_length = distance_along_spine;
                dst_annotations += 2;
            }
            base_index += 2 + spine_length * 2;
        } else {
            dst_annotations -= spine_length * 2;
            for (uint16_t i = 0; i < spine_length; i++) {
                dst_annotations[0].spine_length = distance_along_spine;
                dst_annotations[1].spine_length = distance_along_spine;
                dst_annotations += 2;
            }
            base_index += spine_length * 2;
        }
    }

    assert(src_position - spines.vertices == spines.num_vertices);
    assert(dst_annotations - mesh->vertex_annotations == mesh->num_vertices);
    assert(dst_positions - mesh->vertex_positions == mesh->num_vertices);
    assert(dst_indices - mesh->triangle_indices ==
        mesh->num_triangles * ind_per_tri);

    return mesh;
}

par_streamlines_mesh* par_streamlines_draw_curves_cubic(
    par_streamlines_context* context, par_streamlines_spine_list spines)
{
    return &context->result;
}

par_streamlines_mesh* par_streamlines_draw_curves_quadratic(
    par_streamlines_context* context, par_streamlines_spine_list spines)
{
    return &context->result;
}

par_streamlines_mesh* par_streamlines_draw_streamlines(
    par_streamlines_context* context, par_streamlines_callback func,
    uint32_t frame_index)
{
    return &context->result;
}

#endif // PAR_STREAMLINES_IMPLEMENTATION
#endif // PAR_STREAMLINES_H
