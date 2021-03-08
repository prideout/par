// OCTASPHERE :: https://prideout.net/blog/octasphere
// Tiny malloc-free library that generates triangle meshes for spheres, rounded boxes, and capsules.
//
// This library proffers the following functions:
//
//   - par_octasphere_get_counts
//   - par_octasphere_populate
//
// Usage example:
//
//        /* Specify a 100x100x20 rounded box. */
//        const par_octasphere_config cfg = {
//            .corner_radius = 5,
//            .width = 100,
//            .height = 100,
//            .depth = 20,
//            .num_subdivisions = 3,
//        };
//
//        /* Allocate memory for the mesh and opt-out of normals. */
//        uint32_t num_indices;
//        uint32_t num_vertices;
//        par_octasphere_get_counts(&cfg, &num_indices, &num_vertices);
//        par_octasphere_mesh mesh = {
//            .positions = malloc(sizeof(float) * 3 * num_vertices),
//            .normals = NULL,
//            .texcoords = malloc(sizeof(float) * 2 * num_vertices),
//            .indices = malloc(sizeof(uint16_t) * num_indices),
//        };
//
//        /* Generate vertex coordinates, UV's, and triangle indices. */
//        par_octasphere_populate(&cfg, &mesh);
//
// To generate a sphere: set width, height, and depth to 0 in your configuration.
// To generate a capsule shape: set only two of these dimensions to 0.
//
// Distributed under the MIT License, see bottom of file.

#ifndef PAR_OCTASPHERE_H
#define PAR_OCTASPHERE_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PAR_OCTASPHERE_MAX_SUBDIVISIONS 5

typedef enum {
    PAR_OCTASPHERE_UV_LATLONG,  // Classic sphere mapping with theta-phi.
    PAR_OCTASPHERE_UV_PTEX,     // Each face in the control mesh gets its own patch.
} par_octasphere_uv_mode;

typedef enum {
    PAR_OCTASPHERE_NORMALS_SMOOTH,
} par_octasphere_normals_mode;

typedef struct {
    float corner_radius;
    float width;
    float height;
    float depth;
    int num_subdivisions;
    par_octasphere_uv_mode uv_mode;
    par_octasphere_normals_mode normals_mode;
} par_octasphere_config;

typedef struct {
    float* positions;
    float* normals;
    float* texcoords;
    uint16_t* indices;
    uint32_t num_indices;
    uint32_t num_vertices;
} par_octasphere_mesh;

// Computes the number of indices and vertices for the given octasphere config.
void par_octasphere_get_counts(const par_octasphere_config* config, uint32_t* num_indices,
                               uint32_t* num_vertices);

// Populates a pre-allocated mesh structure with indices and vertices.
void par_octasphere_populate(const par_octasphere_config* config, par_octasphere_mesh* mesh);

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------
#ifdef PAR_OCTASPHERE_IMPLEMENTATION

#include <assert.h>
#include <math.h>
#include <memory.h>  // for memcpy

#define PARO_PI (3.14159265359)
#define PARO_MIN(a, b) (a > b ? b : a)
#define PARO_MAX(a, b) (a > b ? a : b)
#define PARO_CLAMP(v, lo, hi) PARO_MAX(lo, PARO_MIN(hi, v))
#define PARO_MAX_BOUNDARY_LENGTH ((1 << PAR_OCTASPHERE_MAX_SUBDIVISIONS) + 1)

// Writes two triangles into the index buffer according to the given quad corners.
static uint16_t* paro_write_quad(uint16_t* dst, uint16_t a, uint16_t b, uint16_t c, uint16_t d) {
    *dst++ = a;
    *dst++ = b;
    *dst++ = c;
    *dst++ = c;
    *dst++ = d;
    *dst++ = a;
    return dst;
}

// Copies the 4 verts at the given indices, appending them to the end of the vertex buffer.
static void paro_write_quad_unwelded(float const* src_vertices, uint16_t** pindex_write_ptr,
                                     float** pvertex_write_ptr, uint16_t a, uint16_t b, uint16_t c,
                                     uint16_t d) {
    float* vertex_write_ptr = *pvertex_write_ptr;
    uint16_t* index_write_ptr = *pindex_write_ptr;
    const uint16_t first_index = (vertex_write_ptr - src_vertices) / 3;
    *index_write_ptr++ = first_index + 0;
    *index_write_ptr++ = first_index + 1;
    *index_write_ptr++ = first_index + 2;
    *index_write_ptr++ = first_index + 2;
    *index_write_ptr++ = first_index + 3;
    *index_write_ptr++ = first_index + 0;
    *pindex_write_ptr = index_write_ptr;
    for (int axis = 0; axis < 3; ++axis) *vertex_write_ptr++ = src_vertices[a * 3 + axis];
    for (int axis = 0; axis < 3; ++axis) *vertex_write_ptr++ = src_vertices[b * 3 + axis];
    for (int axis = 0; axis < 3; ++axis) *vertex_write_ptr++ = src_vertices[c * 3 + axis];
    for (int axis = 0; axis < 3; ++axis) *vertex_write_ptr++ = src_vertices[d * 3 + axis];
    *pvertex_write_ptr = vertex_write_ptr;
}

static void paro_write_quad_uv(float** puv_write_ptr, float x0, float y0, float x1, float y1, int a,
                               int b, int c, int d) {
    // clang-format off
    float corners[4][2];
    corners[0][0] = x0; corners[0][1] = y0;
    corners[1][0] = x0; corners[1][1] = y1;
    corners[2][0] = x1; corners[2][1] = y0;
    corners[3][0] = x1; corners[3][1] = y1;
    float* uv = *puv_write_ptr;
    *uv++ = corners[a][0]; *uv++ = corners[a][1];
    *uv++ = corners[b][0]; *uv++ = corners[b][1];
    *uv++ = corners[c][0]; *uv++ = corners[c][1];
    *uv++ = corners[c][0]; *uv++ = corners[d][1];
    // clang-format on
    *puv_write_ptr = uv;
}

static void paro_write_ui3(uint16_t* dst, int index, uint16_t a, uint16_t b, uint16_t c) {
    dst[index * 3 + 0] = a;
    dst[index * 3 + 1] = b;
    dst[index * 3 + 2] = c;
}

static float* paro_write_f3(float* dst, const float src[3]) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    return dst + 3;
}

static void paro_copy(float dst[3], const float src[3]) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
}

static float paro_dot(const float a[3], const float b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void paro_add(float result[3], float const a[3], float const b[3]) {
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
}

static void paro_normalize(float v[3]) {
    float lsqr = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (lsqr > 0) {
        v[0] /= lsqr;
        v[1] /= lsqr;
        v[2] /= lsqr;
    }
}

static void paro_cross(float result[3], float const a[3], float const b[3]) {
    float x = (a[1] * b[2]) - (a[2] * b[1]);
    float y = (a[2] * b[0]) - (a[0] * b[2]);
    float z = (a[0] * b[1]) - (a[1] * b[0]);
    result[0] = x;
    result[1] = y;
    result[2] = z;
}

static void paro_scale(float dst[3], float v) {
    dst[0] *= v;
    dst[1] *= v;
    dst[2] *= v;
}

static void paro_scaled(float dst[3], const float src[3], float v) {
    dst[0] = src[0] * v;
    dst[1] = src[1] * v;
    dst[2] = src[2] * v;
}

static void paro_quat_from_rotation(float quat[4], const float axis[3], float radians) {
    paro_copy(quat, axis);
    paro_normalize(quat);
    paro_scale(quat, sin(0.5 * radians));
    quat[3] = cos(0.5 * radians);
}

static void paro_quat_from_eulers(float quat[4], const float eulers[3]) {
    const float roll = eulers[0];
    const float pitch = eulers[1];
    const float yaw = eulers[2];
    const float halfRoll = roll * 0.5;
    const float sR = sin(halfRoll);
    const float cR = cos(halfRoll);
    const float halfPitch = pitch * 0.5;
    const float sP = sin(halfPitch);
    const float cP = cos(halfPitch);
    const float halfYaw = yaw * 0.5;
    const float sY = sin(halfYaw);
    const float cY = cos(halfYaw);
    quat[0] = (sR * cP * cY) + (cR * sP * sY);
    quat[1] = (cR * sP * cY) - (sR * cP * sY);
    quat[2] = (cR * cP * sY) + (sR * sP * cY);
    quat[3] = (cR * cP * cY) - (sR * sP * sY);
}

static void paro_quat_rotate_vector(float dst[3], const float quat[4], const float src[3]) {
    float t[3];
    paro_cross(t, quat, src);
    paro_scale(t, 2.0);

    float p[3];
    paro_cross(p, quat, t);

    paro_scaled(dst, t, quat[3]);
    paro_add(dst, dst, src);
    paro_add(dst, dst, p);
}

static float* paro_write_geodesic(float* dst, const float point_a[3], const float point_b[3],
                                  int num_segments) {
    dst = paro_write_f3(dst, point_a);
    if (num_segments == 0) {
        return dst;
    }
    const float angle_between_endpoints = acos(paro_dot(point_a, point_b));
    const float dtheta = angle_between_endpoints / num_segments;
    float rotation_axis[3], quat[4];
    paro_cross(rotation_axis, point_a, point_b);
    for (int point_index = 1; point_index < num_segments; point_index++, dst += 3) {
        paro_quat_from_rotation(quat, rotation_axis, dtheta * point_index);
        paro_quat_rotate_vector(dst, quat, point_a);
    }
    return paro_write_f3(dst, point_b);
}

// Find the vertex indices along each of the first patch's 3 edges.
static void paro_get_patch_boundaries(const par_octasphere_config* config,
                                      uint16_t boundaries[3][PARO_MAX_BOUNDARY_LENGTH]) {
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    int a = 0, b = 0, c = 0, row;
    uint16_t j0 = 0;
    for (int col_index = 0; col_index < n - 1; col_index++) {
        int col_height = n - 1 - col_index;
        uint16_t j1 = j0 + 1;
        boundaries[0][a++] = j0;
        for (row = 0; row < col_height - 1; row++) {
            if (col_height == n - 1) {
                boundaries[2][c++] = j0 + row;
            }
        }
        if (col_height == n - 1) {
            boundaries[2][c++] = j0 + row;
            boundaries[2][c++] = j1 + row;
        }
        boundaries[1][b++] = j1 + row;
        j0 += col_height + 1;
    }
    boundaries[0][a] = boundaries[1][b] = j0 + row;
}

static void paro_add_quads_ptex(const par_octasphere_config* config, par_octasphere_mesh* mesh) {
    assert(config->uv_mode == PAR_OCTASPHERE_UV_PTEX);
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    const int verts_per_patch = n * (n + 1) / 2;

#warning "TODO: Compute normals for newly added verts."

    // - 4*(n-1) quads between the 4 top patches.
    // - 4*(n-1) quads between the 4 bottom patches.
    // - 4*(n-1) quads between the top and bottom patches.
    // - 6 quads to fill "holes" in each cuboid face.
    const int num_connection_quads = (4 + 4 + 4) * (n - 1) + 6;

    uint16_t boundaries[3][PARO_MAX_BOUNDARY_LENGTH];
    paro_get_patch_boundaries(config, boundaries);

    uint16_t* ind_write_ptr = mesh->indices + mesh->num_indices;
    float* pos_write_ptr = mesh->positions + mesh->num_vertices * 3;

    // Go around the top half.
    for (int patch = 0; patch < 4; patch++) {
        const int next_patch = (patch + 1) % 4;
        const uint16_t* boundary_a = boundaries[1];
        const uint16_t* boundary_b = boundaries[0];
        const uint16_t offset_a = verts_per_patch * patch;
        const uint16_t offset_b = verts_per_patch * next_patch;
        for (int i = 0; i < n - 1; i++) {
            const uint16_t a = boundary_a[i] + offset_a;
            const uint16_t b = boundary_b[i] + offset_b;
            const uint16_t c = boundary_a[i + 1] + offset_a;
            const uint16_t d = boundary_b[i + 1] + offset_b;
            paro_write_quad_unwelded(mesh->positions, &ind_write_ptr, &pos_write_ptr, a, b, d, c);
        }
    }

    // Go around the bottom half.
    for (int patch = 4; patch < 8; patch++) {
        const int next_patch = 4 + (patch + 1) % 4;
        const uint16_t* boundary_a = boundaries[0];
        const uint16_t* boundary_b = boundaries[2];
        const uint16_t offset_a = verts_per_patch * patch;
        const uint16_t offset_b = verts_per_patch * next_patch;
        for (int i = 0; i < n - 1; i++) {
            const uint16_t a = boundary_a[i] + offset_a;
            const uint16_t b = boundary_b[i] + offset_b;
            const uint16_t c = boundary_a[i + 1] + offset_a;
            const uint16_t d = boundary_b[i + 1] + offset_b;
            paro_write_quad_unwelded(mesh->positions, &ind_write_ptr, &pos_write_ptr, d, b, a, c);
        }
    }
    // Connect the top and bottom halves.
    for (int patch = 0; patch < 4; patch++) {
        const int next_patch = 4 + (4 - patch) % 4;
        const uint16_t* boundary_a = boundaries[2];
        const uint16_t* boundary_b = boundaries[1];
        const uint16_t offset_a = verts_per_patch * patch;
        const uint16_t offset_b = verts_per_patch * next_patch;
        for (int i = 0; i < n - 1; i++) {
            const uint16_t a = boundary_a[i] + offset_a;
            const uint16_t b = boundary_b[n - 1 - i] + offset_b;
            const uint16_t c = boundary_a[i + 1] + offset_a;
            const uint16_t d = boundary_b[n - 1 - i - 1] + offset_b;
            paro_write_quad_unwelded(mesh->positions, &ind_write_ptr, &pos_write_ptr, a, b, d, c);
        }
    }

    // Fill in the top and bottom holes.
    uint16_t a, b, c, d;
    a = boundaries[0][n - 1];
    b = a + verts_per_patch;
    c = b + verts_per_patch;
    d = c + verts_per_patch;
    paro_write_quad_unwelded(mesh->positions, &ind_write_ptr, &pos_write_ptr, a, b, c, d);
    a = boundaries[2][0] + verts_per_patch * 4;
    b = a + verts_per_patch;
    c = b + verts_per_patch;
    d = c + verts_per_patch;
    paro_write_quad_unwelded(mesh->positions, &ind_write_ptr, &pos_write_ptr, a, b, c, d);

    // Fill in the side holes.
    const int sides[4][2] = {{7, 0}, {1, 2}, {3, 4}, {5, 6}};
    for (int side = 0; side < 4; side++) {
        int patch_index, patch, next_patch;
        uint16_t *boundary_a, *boundary_b;
        uint16_t offset_a, offset_b;

        uint16_t a, b;
        patch_index = sides[side][0];
        patch = patch_index / 2;
        next_patch = 4 + (4 - patch) % 4;
        offset_a = verts_per_patch * patch;
        offset_b = verts_per_patch * next_patch;
        boundary_a = boundaries[2];
        boundary_b = boundaries[1];
        if (patch_index % 2 == 0) {
            a = boundary_a[0] + offset_a;
            b = boundary_b[n - 1] + offset_b;
        } else {
            a = boundary_a[n - 1] + offset_a;
            b = boundary_b[0] + offset_b;
        }

        uint16_t c, d;
        patch_index = sides[side][1];
        patch = patch_index / 2;
        next_patch = 4 + (4 - patch) % 4;
        offset_a = verts_per_patch * patch;
        offset_b = verts_per_patch * next_patch;
        boundary_a = boundaries[2];
        boundary_b = boundaries[1];
        if (patch_index % 2 == 0) {
            c = boundary_a[0] + offset_a;
            d = boundary_b[n - 1] + offset_b;
        } else {
            c = boundary_a[n - 1] + offset_a;
            d = boundary_b[0] + offset_b;
        }
        paro_write_quad_unwelded(mesh->positions, &ind_write_ptr, &pos_write_ptr, a, b, d, c);
    }

    mesh->num_indices += num_connection_quads * 6;
    mesh->num_vertices += num_connection_quads * 4;

#ifndef NDEBUG
    assert(mesh->num_indices == ind_write_ptr - mesh->indices);
    assert(mesh->num_vertices == (pos_write_ptr - mesh->positions) / 3);
    uint32_t expected_indices;
    uint32_t expected_vertices;
    par_octasphere_get_counts(config, &expected_indices, &expected_vertices);
    assert(mesh->num_indices == expected_indices);
    assert(mesh->num_vertices == expected_vertices);
#endif
}

// This is similar to paro_add_quads_ptex except that verts between patches and quads are welded.
static void paro_add_quads(const par_octasphere_config* config, par_octasphere_mesh* mesh) {
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    const int verts_per_patch = n * (n + 1) / 2;

    uint16_t boundaries[3][PARO_MAX_BOUNDARY_LENGTH];
    paro_get_patch_boundaries(config, boundaries);

    uint16_t* write_ptr = mesh->indices + mesh->num_indices;
    const uint16_t* begin_ptr = write_ptr;

    // Go around the top half.
    for (int patch = 0; patch < 4; patch++) {
        const int next_patch = (patch + 1) % 4;
        const uint16_t* boundary_a = boundaries[1];
        const uint16_t* boundary_b = boundaries[0];
        const uint16_t offset_a = verts_per_patch * patch;
        const uint16_t offset_b = verts_per_patch * next_patch;
        for (int i = 0; i < n - 1; i++) {
            const uint16_t a = boundary_a[i] + offset_a;
            const uint16_t b = boundary_b[i] + offset_b;
            const uint16_t c = boundary_a[i + 1] + offset_a;
            const uint16_t d = boundary_b[i + 1] + offset_b;
            write_ptr = paro_write_quad(write_ptr, a, b, d, c);
        }
    }
    // Go around the bottom half.
    for (int patch = 4; patch < 8; patch++) {
        const int next_patch = 4 + (patch + 1) % 4;
        const uint16_t* boundary_a = boundaries[0];
        const uint16_t* boundary_b = boundaries[2];
        const uint16_t offset_a = verts_per_patch * patch;
        const uint16_t offset_b = verts_per_patch * next_patch;
        for (int i = 0; i < n - 1; i++) {
            const uint16_t a = boundary_a[i] + offset_a;
            const uint16_t b = boundary_b[i] + offset_b;
            const uint16_t c = boundary_a[i + 1] + offset_a;
            const uint16_t d = boundary_b[i + 1] + offset_b;
            write_ptr = paro_write_quad(write_ptr, d, b, a, c);
        }
    }
    // Connect the top and bottom halves.
    for (int patch = 0; patch < 4; patch++) {
        const int next_patch = 4 + (4 - patch) % 4;
        const uint16_t* boundary_a = boundaries[2];
        const uint16_t* boundary_b = boundaries[1];
        const uint16_t offset_a = verts_per_patch * patch;
        const uint16_t offset_b = verts_per_patch * next_patch;
        for (int i = 0; i < n - 1; i++) {
            const uint16_t a = boundary_a[i] + offset_a;
            const uint16_t b = boundary_b[n - 1 - i] + offset_b;
            const uint16_t c = boundary_a[i + 1] + offset_a;
            const uint16_t d = boundary_b[n - 1 - i - 1] + offset_b;
            write_ptr = paro_write_quad(write_ptr, a, b, d, c);
        }
    }

    // Fill in the top and bottom holes.
    uint16_t a, b, c, d;
    a = boundaries[0][n - 1];
    b = a + verts_per_patch;
    c = b + verts_per_patch;
    d = c + verts_per_patch;
    write_ptr = paro_write_quad(write_ptr, a, b, c, d);
    a = boundaries[2][0] + verts_per_patch * 4;
    b = a + verts_per_patch;
    c = b + verts_per_patch;
    d = c + verts_per_patch;
    write_ptr = paro_write_quad(write_ptr, a, b, c, d);

    // Fill in the side holes.
    const int sides[4][2] = {{7, 0}, {1, 2}, {3, 4}, {5, 6}};
    for (int side = 0; side < 4; side++) {
        int patch_index, patch, next_patch;
        uint16_t *boundary_a, *boundary_b;
        uint16_t offset_a, offset_b;

        uint16_t a, b;
        patch_index = sides[side][0];
        patch = patch_index / 2;
        next_patch = 4 + (4 - patch) % 4;
        offset_a = verts_per_patch * patch;
        offset_b = verts_per_patch * next_patch;
        boundary_a = boundaries[2];
        boundary_b = boundaries[1];
        if (patch_index % 2 == 0) {
            a = boundary_a[0] + offset_a;
            b = boundary_b[n - 1] + offset_b;
        } else {
            a = boundary_a[n - 1] + offset_a;
            b = boundary_b[0] + offset_b;
        }

        uint16_t c, d;
        patch_index = sides[side][1];
        patch = patch_index / 2;
        next_patch = 4 + (4 - patch) % 4;
        offset_a = verts_per_patch * patch;
        offset_b = verts_per_patch * next_patch;
        boundary_a = boundaries[2];
        boundary_b = boundaries[1];
        if (patch_index % 2 == 0) {
            c = boundary_a[0] + offset_a;
            d = boundary_b[n - 1] + offset_b;
        } else {
            c = boundary_a[n - 1] + offset_a;
            d = boundary_b[0] + offset_b;
        }
        write_ptr = paro_write_quad(write_ptr, a, b, d, c);
    }

    mesh->num_indices += write_ptr - begin_ptr;

#ifndef NDEBUG
    uint32_t expected_indices;
    uint32_t expected_vertices;
    par_octasphere_get_counts(config, &expected_indices, &expected_vertices);
    assert(mesh->num_indices == expected_indices);
#endif
}

static void paro_populate_uvs_latlong(const par_octasphere_config* config,
                                      par_octasphere_mesh* mesh) {
    assert(config->uv_mode == PAR_OCTASPHERE_UV_LATLONG);
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    const int verts_per_patch = n * (n + 1) / 2;
    const int total_vertices = verts_per_patch * 8;
    for (int i = 0; i < total_vertices; i++) {
        const int octant = i / verts_per_patch;
        const int relative_index = i % verts_per_patch;
        float* uv = mesh->texcoords + i * 2;
        const float* xyz = mesh->positions + i * 3;
        const float x = xyz[0], y = xyz[1], z = xyz[2];
        const float phi = -atan2(z, x);
        const float theta = acos(y);
        uv[0] = 0.5 * (phi / PARO_PI + 1.0);
        uv[1] = theta / PARO_PI;
        // Special case for the north pole.
        if (octant < 4 && relative_index == verts_per_patch - 1) {
            uv[0] = fmod(0.375 + 0.25 * octant, 1.0);
            uv[1] = 0;
        }
        // Special case for the south pole.
        if (octant >= 4 && relative_index == 0) {
            uv[0] = 0.375 - 0.25 * (octant - 4);
            uv[0] = uv[0] + uv[0] < 0 ? 1.0 : 0.0;
            uv[1] = 1.0;
        }
        // Adjust the prime meridian for proper wrapping.
        if ((octant == 2 || octant == 6) && uv[0] < 0.5) {
            uv[0] += 1.0;
        }
    }
}

static void paro_octahedral_proj(const float dir_in[3], float uv_out[2], int octant) {
    float dir[3];
    paro_copy(dir, dir_in);
    paro_scale(dir, 1.0f / (fabs(dir[0]) + fabs(dir[1]) + fabs(dir[2])));
    const float rev[2] = {
        fabs(dir[2]) - 1.0f,
        fabs(dir[0]) - 1.0f,
    };
    const float neg[2] = {
        dir[0] < 0 ? rev[0] : -rev[0],
        dir[2] < 0 ? rev[1] : -rev[1],
    };
    uv_out[0] = dir[1] < 0 ? neg[0] : dir[0];
    uv_out[1] = dir[1] < 0 ? neg[1] : dir[2];
    uv_out[0] = 0.5 * uv_out[0] + 0.5;
    uv_out[1] = 0.5 * uv_out[1] + 0.5;
    if (octant == 4) {  // G
        if (uv_out[1] <= 0.5) uv_out[1] = 1.0 - uv_out[1];
    }
    if (octant == 5) {  // E
        if (uv_out[1] <= 0.5) uv_out[1] = 1.0 - uv_out[1];
        if (uv_out[0] > 0.5) uv_out[0] = 1.0 - uv_out[0];
    }
    if (octant == 6) {  // A
        if (uv_out[0] > 0.5) uv_out[0] = 1.0 - uv_out[0];
    }
    if (octant == 7) {  // C
        if (uv_out[0] <= 0.5) uv_out[0] = 1.0 - uv_out[0];
        if (uv_out[1] > 0.5) uv_out[1] = 1.0 - uv_out[1];
    }
}

static void paro_populate_uvs_ptex_patches(const par_octasphere_config* config,
                                           par_octasphere_mesh* mesh) {
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    const int verts_per_patch = n * (n + 1) / 2;
    assert(config->uv_mode == PAR_OCTASPHERE_UV_PTEX);
    const float* pos_read_ptr = mesh->positions;
    float* uvs_write_ptr = mesh->texcoords;

    for (int i = 0; i < verts_per_patch * 4; i++, pos_read_ptr += 3, uvs_write_ptr += 2) {
        paro_octahedral_proj(pos_read_ptr, uvs_write_ptr, 0);  // B, D, F, H
        uvs_write_ptr[0] *= 2.0 / 11.0;
    }

    for (int i = 0; i < verts_per_patch; i++, pos_read_ptr += 3, uvs_write_ptr += 2) {
        paro_octahedral_proj(pos_read_ptr, uvs_write_ptr, 4);  // G
        uvs_write_ptr[0] *= 2.0 / 11.0;
    }

    for (int i = 0; i < verts_per_patch; i++, pos_read_ptr += 3, uvs_write_ptr += 2) {
        paro_octahedral_proj(pos_read_ptr, uvs_write_ptr, 5);  // E
        uvs_write_ptr[0] *= 2.0 / 11.0;
    }

    for (int i = 0; i < verts_per_patch; i++, pos_read_ptr += 3, uvs_write_ptr += 2) {
        paro_octahedral_proj(pos_read_ptr, uvs_write_ptr, 6);  // A
        uvs_write_ptr[0] *= 2.0 / 11.0;
    }

    for (int i = 0; i < verts_per_patch; i++, pos_read_ptr += 3, uvs_write_ptr += 2) {
        paro_octahedral_proj(pos_read_ptr, uvs_write_ptr, 7);  // C
        uvs_write_ptr[0] *= 2.0 / 11.0;
    }
}

static void paro_populate_uvs_ptex_quads(const par_octasphere_config* config,
                                         par_octasphere_mesh* mesh) {
    assert(config->uv_mode == PAR_OCTASPHERE_UV_PTEX);
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    const int verts_per_patch = n * (n + 1) / 2;
    float* uvs_write_ptr = mesh->texcoords + verts_per_patch * 8 * 2;

    for (int i = verts_per_patch * 8; i < mesh->num_vertices; i++, uvs_write_ptr += 2) {
        uvs_write_ptr[0] = 0;
        uvs_write_ptr[1] = 0;
    }
    uvs_write_ptr = mesh->texcoords + verts_per_patch * 8 * 2;

    float x = 0;
    float y = 0;

    const float dx = 1.0 / (n - 1);
    const int a = 0, b = 1, c = 2, d = 3;

    // Go around the top half. (I K M O)
    for (int patch = 0; patch < 4; patch++) {
        for (int i = 0; i < n - 1; i++) {
            paro_write_quad_uv(&uvs_write_ptr, x, 0, x + dx, 1, a, b, d, c);
            x = x + dx;
        }
    }

    // Go around the bottom half. (Q S U W)
    for (int patch = 0; patch < 4; patch++) {
        for (int i = 0; i < n - 1; i++) {
            paro_write_quad_uv(&uvs_write_ptr, x, y, x + 1, y + dx, d, c, a, b);
            y = y + dx;
        }
        y = 0;
        x++;
    }

    // Connect the top and bottom halves. (Y J L N)
    for (int patch = 0; patch < 4; patch++) {
        for (int i = 0; i < n - 1; i++) {
            paro_write_quad_uv(&uvs_write_ptr, x, 0, x + dx, 1, a, b, d, c);
            x = x + dx;
        }
    }

    // Normalize the coordinates.
    uvs_write_ptr = mesh->texcoords + verts_per_patch * 8 * 2;
    for (int patch = 0; patch < 12; patch++) {
        for (int i = 0; i < 4 * (n - 1); i++) {
            if (patch >= 9) {
                uvs_write_ptr[0] -= 9;
                uvs_write_ptr[1] += 1;
            }
            uvs_write_ptr[0] = 2.0 / 11.0 + uvs_write_ptr[0] / 11.0;
            uvs_write_ptr[1] = uvs_write_ptr[1] / 2.0;
            uvs_write_ptr += 2;
        }
    }

    // Central quads of the six cube faces (P R T V X Z)
    // clang-format off
    paro_write_quad_uv(&uvs_write_ptr, x, 1, x + 1, 2, d, c, a, b); x++;
    paro_write_quad_uv(&uvs_write_ptr, x, 1, x + 1, 2, d, c, a, b); x++;
    paro_write_quad_uv(&uvs_write_ptr, x, 1, x + 1, 2, d, c, a, b); x++;
    paro_write_quad_uv(&uvs_write_ptr, x, 1, x + 1, 2, d, c, a, b); x++;
    paro_write_quad_uv(&uvs_write_ptr, x, 1, x + 1, 2, d, c, a, b); x++;
    paro_write_quad_uv(&uvs_write_ptr, x, 1, x + 1, 2, d, c, a, b); x++;
    // clang-format on

    // Normalize the coordinates.
    uvs_write_ptr -= 4 * 6 * 2;
    for (int i = 0; i < 6 * 4; i++) {
        uvs_write_ptr[0] -= 9;
        uvs_write_ptr[0] = 2.0 / 11.0 + uvs_write_ptr[0] / 11.0;
        uvs_write_ptr[1] = uvs_write_ptr[1] / 2.0;
        uvs_write_ptr += 2;
    }
}

void par_octasphere_get_counts(const par_octasphere_config* config, uint32_t* num_indices,
                               uint32_t* num_vertices) {
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    const int verts_per_patch = n * (n + 1) / 2;
    const int triangles_per_patch = (n - 2) * (n - 1) + n - 1;

    // - 4*(n-1) quads between the 4 top patches.
    // - 4*(n-1) quads between the 4 bottom patches.
    // - 4*(n-1) quads between the top and bottom patches.
    // - 6 quads to fill "holes" in each cuboid face.
    const int num_connection_quads = (4 + 4 + 4) * (n - 1) + 6;

    *num_indices = (triangles_per_patch * 8 + num_connection_quads * 2) * 3;
    *num_vertices = verts_per_patch * 8;

    if (config->uv_mode == PAR_OCTASPHERE_UV_PTEX) {
        *num_vertices += num_connection_quads * 4;
    }
}

void par_octasphere_populate(const par_octasphere_config* config, par_octasphere_mesh* mesh) {
    const int ndivisions = PARO_CLAMP(config->num_subdivisions, 0, PAR_OCTASPHERE_MAX_SUBDIVISIONS);
    const int n = (1 << ndivisions) + 1;
    const int verts_per_patch = n * (n + 1) / 2;
    const float r2 = config->corner_radius * 2;
    const float w = PARO_MAX(config->width, r2);
    const float h = PARO_MAX(config->height, r2);
    const float d = PARO_MAX(config->depth, r2);
    const float tx = (w - r2) / 2, ty = (h - r2) / 2, tz = (d - r2) / 2;
    const int triangles_per_patch = (n - 2) * (n - 1) + n - 1;
    const int total_vertices = verts_per_patch * 8;

    // START TESSELLATION OF SINGLE PATCH (one-eighth of the octasphere)
    float* write_ptr = mesh->positions;
    for (int i = 0; i < n; i++) {
        const float theta = PARO_PI * 0.5 * i / (n - 1);
        const float point_a[] = {0, sinf(theta), cosf(theta)};
        const float point_b[] = {cosf(theta), sinf(theta), 0};
        const int num_segments = n - 1 - i;
        write_ptr = paro_write_geodesic(write_ptr, point_a, point_b, num_segments);
    }
    int f = 0, j0 = 0;
    uint16_t* faces = mesh->indices;
    for (int col_index = 0; col_index < n - 1; col_index++) {
        const int col_height = n - 1 - col_index;
        const int j1 = j0 + 1;
        const int j2 = j0 + col_height + 1;
        const int j3 = j0 + col_height + 2;
        for (int row = 0; row < col_height - 1; row++) {
            paro_write_ui3(faces, f++, j0 + row, j1 + row, j2 + row);
            paro_write_ui3(faces, f++, j2 + row, j1 + row, j3 + row);
        }
        const int row = col_height - 1;
        paro_write_ui3(faces, f++, j0 + row, j1 + row, j2 + row);
        j0 = j2;
    }
    // END TESSELLATION OF SINGLE PATCH

    // START 8-WAY CLONE OF PATCH
    // clang-format off
    float euler_angles[8][3] = {
        {0, 0, 0}, {0, 1, 0}, {0, 2, 0}, {0, 3, 0},
        {1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {1, 0, 3},
    };
    // clang-format on
    for (int octant = 1; octant < 8; octant++) {
        paro_scale(euler_angles[octant], PARO_PI * 0.5);
        float quat[4];
        paro_quat_from_eulers(quat, euler_angles[octant]);
        float* dst = mesh->positions + octant * verts_per_patch * 3;
        const float* src = mesh->positions;
        for (int vindex = 0; vindex < verts_per_patch; vindex++, dst += 3, src += 3) {
            paro_quat_rotate_vector(dst, quat, src);
        }
    }
    for (int octant = 1; octant < 8; octant++) {
        const int indices_per_patch = triangles_per_patch * 3;
        uint16_t* dst = mesh->indices + octant * indices_per_patch;
        const uint16_t* src = mesh->indices;
        const uint16_t offset = verts_per_patch * octant;
        for (int iindex = 0; iindex < indices_per_patch; ++iindex) {
            dst[iindex] = src[iindex] + offset;
        }
    }
    // END 8-WAY CLONE OF PATCH

    if (mesh->texcoords && config->uv_mode == PAR_OCTASPHERE_UV_LATLONG) {
        paro_populate_uvs_latlong(config, mesh);
    }

    // The following memcpy for normals works because these patches are spherical in nature, and we
    // have not yet scaled and translated anything. Note that in PTEX mode, we need to populate
    // additional normals for the quads.
    if (mesh->normals && config->normals_mode == PAR_OCTASPHERE_NORMALS_SMOOTH) {
        memcpy(mesh->normals, mesh->positions, sizeof(float) * 3 * total_vertices);
    }

    if (config->corner_radius != 1.0) {
        for (int i = 0; i < total_vertices; i++) {
            float* xyz = mesh->positions + i * 3;
            xyz[0] *= config->corner_radius;
            xyz[1] *= config->corner_radius;
            xyz[2] *= config->corner_radius;
        }
    }

    mesh->num_indices = triangles_per_patch * 8 * 3;
    mesh->num_vertices = total_vertices;

    if (config->uv_mode == PAR_OCTASPHERE_UV_PTEX && mesh->texcoords) {
        paro_populate_uvs_ptex_patches(config, mesh);
    }

    for (int i = 0; i < total_vertices; i++) {
        float* xyz = mesh->positions + i * 3;
        const int octant = i / verts_per_patch;
        const float sx = (octant < 2 || octant == 4 || octant == 7) ? +1 : -1;
        const float sy = octant < 4 ? +1 : -1;
        const float sz = (octant == 0 || octant == 3 || octant == 4 || octant == 5) ? +1 : -1;
        xyz[0] += tx * sx;
        xyz[1] += ty * sy;
        xyz[2] += tz * sz;
    }

    if (config->uv_mode == PAR_OCTASPHERE_UV_PTEX) {
        paro_add_quads_ptex(config, mesh);
        if (mesh->texcoords) {
            paro_populate_uvs_ptex_quads(config, mesh);
        }
    } else {
        paro_add_quads(config, mesh);
    }
}

#endif  // PAR_OCTASPHERE_IMPLEMENTATION
#endif  // PAR_OCTASPHERE_H

// par_octasphere is distributed under the MIT license:
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
