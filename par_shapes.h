// SHAPES :: https:github.com/prideout/par
// Mesh generator for parametric surfaces and other simple geometry.
//
//     http://github.prideout.net/c-shapes/
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

#include <stdint.h>
#include <stdbool.h>

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

#define PAR_SHAPES_VERSION 0.0.0

typedef struct par_shapes_mesh_s {
    float* points;
    int npoints;
    uint16_t* triangles;
    int ntriangles;
    float* normals;
    float* tcoords;
} par_shapes_mesh;

#define PAR_SHAPES_SMOOTH_NORMALS (1 << 0)
#define PAR_SHAPES_TEXTURE_COORDS (1 << 2)

char const * const * par_shapes_list_parametric();
par_shapes_mesh* par_shapes_create_parametric(char const*, int slices,
    int stacks, int flags);
par_shapes_mesh* par_shapes_create_disk(float radius, int slices,
    float const* center, float const* normal, int flags);
void par_shapes_free(par_shapes_mesh*);
void par_shapes_export(par_shapes_mesh const*, char const* objfile);
void par_shapes_merge(par_shapes_mesh* dst, par_shapes_mesh const* src);
void par_shapes_translate(par_shapes_mesh*, float x, float y, float z);
void par_shapes_rotate(par_shapes_mesh*, float radians, float const* axis);
void par_shapes_scale(par_shapes_mesh*, float x, float y, float z);

// Reverse the winding of a run of faces.  Useful when drawing the inside of
// a Cornell Box.  Pass 0 for nfaces to reverse every face in the mesh.
void par_shapes_invert(par_shapes_mesh*, int startface, int nfaces);

// Dereference the entire index buffer and replace the point list.
// This creates an inefficient structure, but is useful for drawing facets.
void par_shapes_unweld(par_shapes_mesh* mesh, bool create_indices);

// Consume an unwelded mesh and insert facet normals into the mesh.
void par_shapes_compute_facet_normals(par_shapes_mesh* m);

// Generate points for a 20-sided polyhedron that fits in the unit sphere.
par_shapes_mesh* par_shapes_create_icosahedron();

// Generate points for a 12-sided polyhedron that fits in the unit sphere.
par_shapes_mesh* par_shapes_create_dodecahedron();

// Create a sphere from a subdivided icosahedron without normals or uvs.
par_shapes_mesh* par_shapes_create_sphere(int nsubdivisions);

// Create a rock shape that sits on the Y=0 plane, and sinks into it a bit.
par_shapes_mesh* par_shapes_create_rock(int seed, int nsubdivisions);

// Create a crappy cloud shape that floats in the Y=0 plane.
par_shapes_mesh* par_shapes_create_cloud(int seed, int nsubdivisions);

// TBD, http://prideout.net/blog/?p=44
typedef void (*par_shapes_fn)(float* const, float*);
par_shapes_mesh* par_shapes_create_custom_parametric(par_shapes_fn, int slices,
    int stacks, int flags);
par_shapes_mesh* par_shapes_create_tree(int seed, int flags);
par_shapes_mesh* par_shapes_create_octohedron();
par_shapes_mesh* par_shapes_create_cube(); // for Cornell boxes

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_SHAPES_IMPLEMENTATION

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_MALLOC(T, N) ((T*) malloc(N * sizeof(T)))
#define PAR_CALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#define PAR_PI (3.14159265359)

static void par_shapes_private_sphere(float* const uv, float* xyz);
static void par_shapes_private_plane(float* const uv, float* xyz);
static void par_shapes_private_klein(float* const uv, float* xyz);
static void par_shapes_private_cylinder(float* const uv, float* xyz);
static void par_shapes_private_torus(float* const uv, float* xyz);

struct osn_context;
static int par_simplex_noise(int64_t seed, struct osn_context** ctx);
static void par_simplex_noise_free(struct osn_context* ctx);
static double par_simplex_noise2(struct osn_context* ctx, double x, double y);

static par_shapes_fn par_shapes_functions[] = {
    par_shapes_private_sphere,
    par_shapes_private_plane,
    par_shapes_private_klein,
    par_shapes_private_cylinder,
    par_shapes_private_torus,
};

static const char* par_shapes_names[] = {
    "sphere",
    "plane",
    "klein",
    "cylinder",
    "torus",
};

static void par_shapes_cross3(float* result, float const* a, float const* b)
{
    float x = (a[1] * b[2]) - (a[2] * b[1]);
    float y = (a[2] * b[0]) - (a[0] * b[2]);
    float z = (a[0] * b[1]) - (a[1] * b[0]);
    result[0] = x;
    result[1] = y;
    result[2] = z;
}

static void par_shapes_mix3(float* dst, float const* a, float const* b, float t)
{
    float x = b[0] * t + a[0] * (1 - t);
    float y = b[1] * t + a[1] * (1 - t);
    float z = b[2] * t + a[2] * (1 - t);
    dst[0] = x;
    dst[1] = y;
    dst[2] = z;
}

static void par_shapes_scale3(float* result, float a)
{
    result[0] *= a;
    result[1] *= a;
    result[2] *= a;
}

static void par_shapes_normalize3(float* v)
{
    float lsqr = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (lsqr > 0) {
        par_shapes_scale3(v, 1.0f / lsqr);
    }
}

static void par_shapes_subtract3(float* result, float const* a)
{
    result[0] -= a[0];
    result[1] -= a[1];
    result[2] -= a[2];
}

static void par_shapes_add3(float* result, float const* a)
{
    result[0] += a[0];
    result[1] += a[1];
    result[2] += a[2];
}

char const * const * par_shapes_list_parametric()
{
    return par_shapes_names;
}

par_shapes_mesh* par_shapes_create_parametric(char const* name,
    int slices, int stacks, int flags)
{
    if (slices < 3 || stacks < 3) {
        return 0;
    }
    char const * const * list = par_shapes_names;
    int shape_index = 0;
    while (*list) {
        if (!strcmp(*list, name)) {
            break;
        }
        shape_index++;
        list++;
    }
    if (!*list) {
        return 0;
    }
    par_shapes_fn fn = par_shapes_functions[shape_index];
    par_shapes_mesh* mesh = (par_shapes_mesh*)
        calloc(sizeof(par_shapes_mesh), 1);

    // Generate verts.
    mesh->npoints = (slices + 1) * (stacks + 1);
    mesh->points = PAR_CALLOC(float, 3 * mesh->npoints);
    float uv[2];
    float xyz[3];
    float* points = mesh->points;
    for (int slice = 0; slice < slices + 1; slice++) {
        uv[1] = (float) slice / slices;
        for (int stack = 0; stack < stacks + 1; stack++) {
            uv[0] = (float) stack / stacks;
            fn(uv, xyz);
            *points++ = xyz[0];
            *points++ = xyz[1];
            *points++ = xyz[2];
        }
    }

    // Generate smooth normals.
    if (flags & PAR_SHAPES_SMOOTH_NORMALS) {
        mesh->normals = PAR_CALLOC(float, 3 * mesh->npoints);
        float* normals = mesh->normals;
        float du[2];
        float dv[2];
        float du_xyz[3];
        float dv_xyz[3];
        float epsilon = 0.25 / PAR_MAX(slices, stacks);
        for (int slice = 0; slice < slices + 1; slice++) {
            du[1] = dv[1] = uv[1] = (float) slice / slices;
            for (int stack = 0; stack < stacks + 1; stack++) {
                du[0] = dv[0] = uv[0] = (float) stack / stacks;
                du[0] += epsilon;
                dv[1] += epsilon;
                fn(uv, xyz);
                fn(du, du_xyz);
                fn(dv, dv_xyz);
                par_shapes_subtract3(du_xyz, xyz);
                par_shapes_subtract3(dv_xyz, xyz);
                par_shapes_cross3(normals, du_xyz, dv_xyz);
                par_shapes_normalize3(normals);
                normals += 3;
            }
        }
    }

    // Generate texture coordinates.
    if (flags & PAR_SHAPES_TEXTURE_COORDS) {
        mesh->tcoords = PAR_CALLOC(float, 2 * mesh->npoints);
        float* uvs = mesh->tcoords;
        for (int slice = 0; slice < slices + 1; slice++) {
            uv[1] = (float) slice / slices;
            for (int stack = 0; stack < stacks + 1; stack++) {
                uv[0] = (float) stack / stacks;
                *uvs++ = uv[0];
                *uvs++ = uv[1];
            }
        }
    }

    // Generate faces.
    mesh->ntriangles = 2 * slices * stacks;
    mesh->triangles = (uint16_t*)
        calloc(sizeof(uint16_t) * 3 * mesh->ntriangles, 1);
    int v = 0;
    uint16_t* face = mesh->triangles;
    for (int slice = 0; slice < slices; slice++) {
        for (int stack = 0; stack < stacks; stack++) {
            int next = stack + 1;
            *face++ = v + stack;
            *face++ = v + next;
            *face++ = v + stack + stacks + 1;
            *face++ = v + next;
            *face++ = v + next + stacks + 1;
            *face++ = v + stack + stacks + 1;
        }
        v += stacks + 1;
    }

    // Hacks for single-sided surfaces go here.  :)
    if (!strcmp(name, "klein")) {
        int face = 0;
        for (int slice = 0; slice < slices; slice++) {
            for (int stack = 0; stack < stacks; stack++, face += 2) {
                if (stack < 27 * stacks / 32) {
                    par_shapes_invert(mesh, face, 2);
                }
            }
        }
    }

    return mesh;
}

void par_shapes_free(par_shapes_mesh* mesh)
{
    free(mesh->points);
    free(mesh->triangles);
    free(mesh->normals);
    free(mesh->tcoords);
    free(mesh);
}

void par_shapes_export(par_shapes_mesh const* mesh, char const* filename)
{
    FILE* objfile = fopen(filename, "wt");
    float const* points = mesh->points;
    float const* tcoords = mesh->tcoords;
    float const* norms = mesh->normals;
    uint16_t const* indices = mesh->triangles;
    if (tcoords && norms) {
        for (int nvert = 0; nvert < mesh->npoints; nvert++) {
            fprintf(objfile, "v %f %f %f\n", points[0], points[1], points[2]);
            fprintf(objfile, "vt %f %f\n", tcoords[0], tcoords[1]);
            fprintf(objfile, "vn %f %f %f\n", norms[0], norms[1], norms[2]);
            points += 3;
            norms += 3;
            tcoords += 2;
        }
        for (int nface = 0; nface < mesh->ntriangles; nface++) {
            int a = 1 + *indices++;
            int b = 1 + *indices++;
            int c = 1 + *indices++;
            fprintf(objfile, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                a, a, a, b, b, b, c, c, c);
        }
    } else if (norms) {
        for (int nvert = 0; nvert < mesh->npoints; nvert++) {
            fprintf(objfile, "v %f %f %f\n", points[0], points[1], points[2]);
            fprintf(objfile, "vn %f %f %f\n", norms[0], norms[1], norms[2]);
            points += 3;
            norms += 3;
        }
        for (int nface = 0; nface < mesh->ntriangles; nface++) {
            int a = 1 + *indices++;
            int b = 1 + *indices++;
            int c = 1 + *indices++;
            fprintf(objfile, "f %d//%d %d//%d %d//%d\n", a, a, b, b, c, c);
        }
    } else if (tcoords) {
        for (int nvert = 0; nvert < mesh->npoints; nvert++) {
            fprintf(objfile, "v %f %f %f\n", points[0], points[1], points[2]);
            fprintf(objfile, "vt %f %f\n", tcoords[0], tcoords[1]);
            points += 3;
            tcoords += 2;
        }
        for (int nface = 0; nface < mesh->ntriangles; nface++) {
            int a = 1 + *indices++;
            int b = 1 + *indices++;
            int c = 1 + *indices++;
            fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
    } else {
        for (int nvert = 0; nvert < mesh->npoints; nvert++) {
            fprintf(objfile, "v %f %f %f\n", points[0], points[1], points[2]);
            points += 3;
        }
        for (int nface = 0; nface < mesh->ntriangles; nface++) {
            int a = 1 + *indices++;
            int b = 1 + *indices++;
            int c = 1 + *indices++;
            fprintf(objfile, "f %d %d %d\n", a, b, c);
        }
    }
    fclose(objfile);
}

static void par_shapes_private_sphere(float* const uv, float* xyz)
{
    float phi = uv[0] * PAR_PI;
    float theta = uv[1] * 2 * PAR_PI;
    xyz[0] = cosf(theta) * sinf(phi);
    xyz[1] = sinf(theta) * sinf(phi);
    xyz[2] = cosf(phi);
}

static void par_shapes_private_plane(float* const uv, float* xyz)
{
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = 0;
}

static void par_shapes_private_klein(float* const uv, float* xyz)
{
    float u = uv[0] * PAR_PI;
    float v = uv[1] * 2 * PAR_PI;
    u = u * 2;
    if (u < PAR_PI) {
        xyz[0] = 3 * cosf(u) * (1 + sinf(u)) + (2 * (1 - cosf(u) / 2)) *
            cosf(u) * cosf(v);
        xyz[2] = -8 * sinf(u) - 2 * (1 - cosf(u) / 2) * sinf(u) * cosf(v);
    } else {
        xyz[0] = 3 * cosf(u) * (1 + sinf(u)) + (2 * (1 - cosf(u) / 2)) *
            cosf(v + PAR_PI);
        xyz[2] = -8 * sinf(u);
    }
    xyz[1] = -2 * (1 - cosf(u) / 2) * sinf(v);
}

static void par_shapes_private_cylinder(float* const uv, float* xyz)
{
    float theta = uv[1] * 2 * PAR_PI;
    xyz[0] = sinf(theta);
    xyz[1] = cosf(theta);
    xyz[2] = uv[0];
}

static void par_shapes_private_torus(float* const uv, float* xyz)
{
    float theta = uv[0] * 2 * PAR_PI;
    float phi = uv[1] * 2 * PAR_PI;
    const float major = 1;
    const float minor = 0.2;
    float beta = major + minor * cosf(phi);
    xyz[0] = cosf(theta) * beta;
    xyz[1] = sinf(theta) * beta;
    xyz[2] = sinf(phi) * minor;
}

void par_shapes_merge(par_shapes_mesh* dst, par_shapes_mesh const* src)
{
    uint16_t offset = dst->npoints;
    int npoints = dst->npoints + src->npoints;
    int vecsize = sizeof(float) * 3;
    float* points = PAR_MALLOC(float, vecsize * npoints / 4);
    memcpy(points, dst->points, vecsize * dst->npoints);
    memcpy(points + 3 * dst->npoints, src->points, vecsize * src->npoints);
    free(dst->points);
    dst->points = points;
    dst->npoints = npoints;
    if (src->normals || dst->normals) {
        float* norms = PAR_CALLOC(float, vecsize * npoints / 4);
        if (dst->normals) {
            memcpy(norms, dst->normals, vecsize * offset);
        }
        if (src->normals) {
            memcpy(norms + 3 * offset, src->normals, vecsize * src->npoints);
        }
        free(dst->normals);
        dst->normals = norms;
    }
    if (src->tcoords || dst->tcoords) {
        int uvsize = sizeof(float) * 2;
        float* uvs = PAR_CALLOC(float, uvsize * npoints /4);
        if (dst->tcoords) {
            memcpy(uvs, dst->tcoords, uvsize * offset);
        }
        if (src->tcoords) {
            memcpy(uvs + 2 * offset, src->tcoords, uvsize * src->npoints);
        }
        free(dst->tcoords);
        dst->tcoords = uvs;
    }
    int ntriangles = dst->ntriangles + src->ntriangles;
    int trisize = sizeof(uint16_t) * 3;
    uint16_t* triangles = PAR_MALLOC(uint16_t, trisize * ntriangles / 2);
    memcpy(triangles, dst->triangles, trisize * dst->ntriangles);
    uint16_t* ptriangles = triangles + 3 * dst->ntriangles;
    uint16_t const* striangles = src->triangles;
    for (int i = 0; i < src->ntriangles; i++) {
        *ptriangles++ = offset + *striangles++;
        *ptriangles++ = offset + *striangles++;
        *ptriangles++ = offset + *striangles++;
    }
    free(dst->triangles);
    dst->triangles = triangles;
    dst->ntriangles = ntriangles;
}

par_shapes_mesh* par_shapes_create_disk(float radius, int slices,
    float const* center, float const* normal, int flags)
{
    if (flags & PAR_SHAPES_TEXTURE_COORDS) {
        // No texture coordinates since this doesn't have a 2D domain.
        return 0;
    }
    par_shapes_mesh* mesh = (par_shapes_mesh*)
        calloc(sizeof(par_shapes_mesh), 1);
    mesh->npoints = slices + 1;
    mesh->points = PAR_MALLOC(float, 3 * mesh->npoints);
    float* points = mesh->points;
    *points++ = 0;
    *points++ = 0;
    *points++ = 0;
    for (int i = 0; i < slices; i++) {
        float theta = i * PAR_PI * 2 / slices;
        *points++ = radius * cos(theta);
        *points++ = radius * sin(theta);
        *points++ = 0;
    }
    float nnormal[3] = {normal[0], normal[1], normal[2]};
    par_shapes_normalize3(nnormal);
    if (flags & PAR_SHAPES_SMOOTH_NORMALS) {
        mesh->normals = PAR_MALLOC(float, 3 * mesh->npoints);
        float* norms = mesh->normals;
        for (int i = 0; i < mesh->npoints; i++) {
            *norms++ = nnormal[0];
            *norms++ = nnormal[1];
            *norms++ = nnormal[2];
        }
    }
    mesh->ntriangles = slices;
    mesh->triangles = (uint16_t*)
        malloc(sizeof(uint16_t) * 3 * mesh->ntriangles);
    uint16_t* triangles = mesh->triangles;
    for (int i = 0; i < slices; i++) {
        *triangles++ = 0;
        *triangles++ = 1 + i;
        *triangles++ = 1 + (i + 1) % slices;
    }
    float k[3] = {0, 0, -1};
    float axis[3];
    par_shapes_cross3(axis, nnormal, k);
    par_shapes_normalize3(axis);
    par_shapes_rotate(mesh, acos(nnormal[2]), axis);
    par_shapes_translate(mesh, center[0], center[1], center[2]);
    return mesh;
}

void par_shapes_translate(par_shapes_mesh* m, float x, float y, float z)
{
    float* points = m->points;
    for (int i = 0; i < m->npoints; i++) {
        *points++ += x;
        *points++ += y;
        *points++ += z;
    }
}

void par_shapes_rotate(par_shapes_mesh* mesh, float radians, float const* axis)
{
    float s = sinf(radians);
    float c = cosf(radians);
    float x = axis[0];
    float y = axis[1];
    float z = axis[2];
    float xy = x * y;
    float yz = y * z;
    float zx = z * x;
    float oneMinusC = 1.0f - c;
    float col0[3] = {(((x * x) * oneMinusC) + c),
        ((xy * oneMinusC) + (z * s)), ((zx * oneMinusC) - (y * s))};
    float col1[3] = {((xy * oneMinusC) - (z * s)),
        (((y * y) * oneMinusC) + c), ((yz * oneMinusC) + (x * s))};
    float col2[3] = {((zx * oneMinusC) + (y * s)),
        ((yz * oneMinusC) - (x * s)), (((z * z) * oneMinusC) + c)};
    float* p = mesh->points;
    for (int i = 0; i < mesh->npoints; i++) {
        float x = col0[0] * p[0] + col1[0] * p[1] + col2[0] * p[2];
        float y = col0[1] * p[0] + col1[1] * p[1] + col2[1] * p[2];
        float z = col0[2] * p[0] + col1[2] * p[1] + col2[2] * p[2];
        *p++ = x;
        *p++ = y;
        *p++ = z;
    }
}

void par_shapes_scale(par_shapes_mesh* m, float x, float y, float z)
{
    float* points = m->points;
    for (int i = 0; i < m->npoints; i++) {
        *points++ *= x;
        *points++ *= y;
        *points++ *= z;
    }
}

void par_shapes_invert(par_shapes_mesh* m, int face, int nfaces)
{
    nfaces = nfaces ? nfaces : m->ntriangles;
    uint16_t* tri = m->triangles + face * 3;
    for (int i = 0; i < nfaces; i++) {
        PAR_SWAP(uint16_t, tri[0], tri[2]);
        tri += 3;
    }
}

par_shapes_mesh* par_shapes_create_icosahedron()
{
    static float verts[] = {
         0.000,  0.000,  1.000,
         0.894,  0.000,  0.447,
         0.276,  0.851,  0.447,
        -0.724,  0.526,  0.447,
        -0.724, -0.526,  0.447,
         0.276, -0.851,  0.447,
         0.724,  0.526, -0.447,
        -0.276,  0.851, -0.447,
        -0.894,  0.000, -0.447,
        -0.276, -0.851, -0.447,
         0.724, -0.526, -0.447,
         0.000,  0.000, -1.000
    };
    static uint16_t faces[] = {
        0,1,2,
        0,2,3,
        0,3,4,
        0,4,5,
        0,5,1,
        7,6,11,
        8,7,11,
        9,8,11,
        10,9,11,
        6,10,11,
        6,2,1,
        7,3,2,
        8,4,3,
        9,5,4,
        10,1,5,
        6,7,2,
        7,8,3,
        8,9,4,
        9,10,5,
        10,6,1
    };
    par_shapes_mesh* mesh = (par_shapes_mesh*)
        calloc(sizeof(par_shapes_mesh), 1);
    mesh->npoints = sizeof(verts) / sizeof(verts[0]) / 3;
    mesh->points = PAR_MALLOC(float, sizeof(verts) / 4);
    memcpy(mesh->points, verts, sizeof(verts));
    mesh->ntriangles = sizeof(faces) / sizeof(faces[0]) / 3;
    mesh->triangles = PAR_MALLOC(uint16_t, sizeof(faces) / 2);
    memcpy(mesh->triangles, faces, sizeof(faces));
    return mesh;
}

par_shapes_mesh* par_shapes_create_dodecahedron()
{
    static float verts[20 * 3] = {
        0.607, 0.000, 0.795,
        0.188, 0.577, 0.795,
        -0.491, 0.357, 0.795,
        -0.491, -0.357, 0.795,
        0.188, -0.577, 0.795,
        0.982, 0.000, 0.188,
        0.304, 0.934, 0.188,
        -0.795, 0.577, 0.188,
        -0.795, -0.577, 0.188,
        0.304, -0.934, 0.188,
        0.795, 0.577, -0.188,
        -0.304, 0.934, -0.188,
        -0.982, 0.000, -0.188,
        -0.304, -0.934, -0.188,
        0.795, -0.577, -0.188,
        0.491, 0.357, -0.795,
        -0.188, 0.577, -0.795,
        -0.607, 0.000, -0.795,
        -0.188, -0.577, -0.795,
        0.491, -0.357, -0.795,
    };
    static uint16_t pentagons[12 * 5] = {
        0,1,2,3,4,
        5,10,6,1,0,
        6,11,7,2,1,
        7,12,8,3,2,
        8,13,9,4,3,
        9,14,5,0,4,
        15,16,11,6,10,
        16,17,12,7,11,
        17,18,13,8,12,
        18,19,14,9,13,
        19,15,10,5,14,
        19,18,17,16,15
    };
    int npentagons = sizeof(pentagons) / sizeof(pentagons[0]) / 5;
    par_shapes_mesh* mesh = PAR_CALLOC(par_shapes_mesh, 1);
    int ncorners = sizeof(verts) / sizeof(verts[0]) / 3;
    mesh->npoints = ncorners;
    mesh->points = PAR_MALLOC(float, mesh->npoints * 3);
    memcpy(mesh->points, verts, sizeof(verts));
    uint16_t const* pentagon = pentagons;
    mesh->ntriangles = npentagons * 3;
    mesh->triangles = PAR_MALLOC(uint16_t, mesh->ntriangles * 3);
    uint16_t* tris = mesh->triangles;
    for (int p = 0; p < npentagons; p++, pentagon += 5) {
        *tris++ = pentagon[0];
        *tris++ = pentagon[1];
        *tris++ = pentagon[2];
        *tris++ = pentagon[0];
        *tris++ = pentagon[2];
        *tris++ = pentagon[3];
        *tris++ = pentagon[0];
        *tris++ = pentagon[3];
        *tris++ = pentagon[4];
    }
    return mesh;
}

void par_shapes_unweld(par_shapes_mesh* mesh, bool create_indices)
{
    int npoints = mesh->ntriangles * 3;
    float* points = PAR_MALLOC(float, 3 * npoints);
    float* dst = points;
    uint16_t const* index = mesh->triangles;
    for (int i = 0; i < npoints; i++) {
        float const* src = mesh->points + 3 * (*index++);
        *dst++ = src[0];
        *dst++ = src[1];
        *dst++ = src[2];
    }
    free(mesh->points);
    mesh->points = points;
    mesh->npoints = npoints;
    if (create_indices) {
        uint16_t* triangles = (uint16_t*)
            malloc(sizeof(uint16_t) * 3 * mesh->ntriangles);
        uint16_t* index = triangles;
        for (int i = 0; i < mesh->ntriangles * 3; i++) {
            *index++ = i;
        }
        free(mesh->triangles);
        mesh->triangles = triangles;
    }
}

void par_shapes_compute_facet_normals(par_shapes_mesh* mesh)
{
    assert(mesh->npoints == mesh->ntriangles * 3 && "Must be unwelded.");
    if (mesh->normals) {
        free(mesh->normals);
    }
    mesh->normals = PAR_MALLOC(float, 3 * mesh->npoints);
    float const* p = mesh->points;
    float* n = mesh->normals;
    for (int t = 0; t < mesh->ntriangles; t++, p += 9, n += 9) {
        n[0] = p[0];
        n[1] = p[1];
        n[2] = p[2];
        n[3] = p[3];
        n[4] = p[4];
        n[5] = p[5];
        par_shapes_subtract3(n, p + 6);
        par_shapes_subtract3(n + 3, p + 6);
        par_shapes_cross3(n, n, n + 3);
        par_shapes_normalize3(n);
        n[3] = n[6] = n[0];
        n[4] = n[7] = n[1];
        n[5] = n[8] = n[2];
    }
}

static void par_shapes_subdivide(par_shapes_mesh* mesh)
{
    assert(mesh->npoints == mesh->ntriangles * 3 && "Must be unwelded.");
    int ntriangles = mesh->ntriangles * 4;
    int npoints = ntriangles * 3;
    float* points = PAR_CALLOC(float, npoints * 3);
    float* dpoint = points;
    float const* spoint = mesh->points;
    for (int t = 0; t < mesh->ntriangles; t++, spoint += 9, dpoint += 3) {
        float const* a = spoint;
        float const* b = spoint + 3;
        float const* c = spoint + 6;
        float const* p0 = dpoint;
        float const* p1 = dpoint + 3;
        float const* p2 = dpoint + 6;
        par_shapes_mix3(dpoint, a, b, 0.5);
        par_shapes_mix3(dpoint += 3, b, c, 0.5);
        par_shapes_mix3(dpoint += 3, a, c, 0.5);
        par_shapes_add3(dpoint += 3, a);
        par_shapes_add3(dpoint += 3, p0);
        par_shapes_add3(dpoint += 3, p2);
        par_shapes_add3(dpoint += 3, p0);
        par_shapes_add3(dpoint += 3, b);
        par_shapes_add3(dpoint += 3, p1);
        par_shapes_add3(dpoint += 3, p2);
        par_shapes_add3(dpoint += 3, p1);
        par_shapes_add3(dpoint += 3, c);
    }
    free(mesh->points);
    mesh->points = points;
    mesh->npoints = npoints;
    mesh->ntriangles = ntriangles;
}

par_shapes_mesh* par_shapes_create_sphere(int nsubd)
{
    par_shapes_mesh* mesh = par_shapes_create_icosahedron();
    par_shapes_unweld(mesh, false);
    free(mesh->triangles);
    mesh->triangles = 0;
    while (nsubd--) {
        par_shapes_subdivide(mesh);
    }
    for (int i = 0; i < mesh->npoints; i++) {
        par_shapes_normalize3(mesh->points + i * 3);
    }
    mesh->triangles = PAR_MALLOC(uint16_t, 3 * mesh->ntriangles);
    for (int i = 0; i < mesh->ntriangles * 3; i++) {
        mesh->triangles[i] = i;
    }
    return mesh;
}

par_shapes_mesh* par_shapes_create_rock(int seed, int subd)
{
    par_shapes_mesh* mesh = par_shapes_create_sphere(subd);
    struct osn_context* ctx;
    par_simplex_noise(seed, &ctx);
    for (int p = 0; p < mesh->npoints; p++) {
        float* pt = mesh->points + p * 3;
        float a = 0.25, f = 1.0;
        double n = a * par_simplex_noise2(ctx, f * pt[0], f * pt[2]);
        a *= 0.5; f *= 2;
        n += a * par_simplex_noise2(ctx, f * pt[0], f * pt[2]);
        par_shapes_scale3(pt, 1.0 + n);
        if (pt[1] < 0) {
            pt[1] = -pow(-pt[1], 0.5) / 2;
        }
    }
    par_simplex_noise_free(ctx);
    return mesh;
}

// This is crap.  It should probably use Worley noise or something.  Want to
// improve it? Make a pull request!  For inspiration, see "The Real Time
// Volumetric Cloudscapes of Horizon Zero Dawn".
par_shapes_mesh* par_shapes_create_cloud(int seed, int nsubd)
{
    par_shapes_mesh* mesh = par_shapes_create_icosahedron();
    par_shapes_unweld(mesh, false);
    free(mesh->triangles);
    mesh->triangles = 0;
    while (nsubd--) {
        par_shapes_subdivide(mesh);
    }
    for (int i = 0; i < mesh->npoints; i++) {
        float* v = mesh->points + i * 3;
        float lsqr = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (lsqr > 0) {
            par_shapes_scale3(v, 1.0f / lsqr);
        }
    }
    mesh->triangles = PAR_MALLOC(uint16_t, 3 * mesh->ntriangles);
    for (int i = 0; i < mesh->ntriangles * 3; i++) {
        mesh->triangles[i] = i;
    }
    struct osn_context* ctx;
    par_simplex_noise(seed, &ctx);
    for (int p = 0; p < mesh->npoints; p++) {
        float* pt = mesh->points + p * 3;
        float a = 0.25, f = 1.0;
        double n = a * par_simplex_noise2(ctx, f * pt[0], f * pt[2] + pt[1]);
        a *= 0.5; f *= 2;
        n += a * par_simplex_noise2(ctx, f * pt[0], f * pt[2] + pt[1]);
        pt[1] /= 3.0;
        par_shapes_scale3(pt, 1.0 + n);
    }
    par_simplex_noise_free(ctx);
    return mesh;
}

// -----------------------------------------------------------------------------
// BEGIN OPEN SIMPLEX NOISE
// -----------------------------------------------------------------------------

#define STRETCH_CONSTANT_2D (-0.211324865405187)  // (1 / sqrt(2 + 1) - 1 ) / 2;
#define SQUISH_CONSTANT_2D (0.366025403784439)  // (sqrt(2 + 1) -1) / 2;
#define STRETCH_CONSTANT_3D (-1.0 / 6.0)  // (1 / sqrt(3 + 1) - 1) / 3;
#define SQUISH_CONSTANT_3D (1.0 / 3.0)  // (sqrt(3+1)-1)/3;
#define STRETCH_CONSTANT_4D (-0.138196601125011)  // (1 / sqrt(4 + 1) - 1) / 4;
#define SQUISH_CONSTANT_4D (0.309016994374947)  // (sqrt(4 + 1) - 1) / 4;

#define NORM_CONSTANT_2D (47.0)
#define NORM_CONSTANT_3D (103.0)
#define NORM_CONSTANT_4D (30.0)

#define DEFAULT_SEED (0LL)

struct osn_context {
    int16_t* perm;
    int16_t* permGradIndex3D;
};

#define ARRAYSIZE(x) (sizeof((x)) / sizeof((x)[0]))

/*
 * Gradients for 2D. They approximate the directions to the
 * vertices of an octagon from the center.
 */
static const int8_t gradients2D[] = {
    5, 2, 2, 5, -5, 2, -2, 5, 5, -2, 2, -5, -5, -2, -2, -5,
};

/*
 * Gradients for 3D. They approximate the directions to the
 * vertices of a rhombicuboctahedron from the center, skewed so
 * that the triangular and square facets can be inscribed inside
 * circles of the same radius.
 */
static const signed char gradients3D[] = {
    -11, 4, 4, -4, 11, 4, -4, 4, 11, 11, 4, 4, 4, 11, 4, 4, 4, 11, -11, -4, 4,
    -4, -11, 4, -4, -4, 11, 11, -4, 4, 4, -11, 4, 4, -4, 11, -11, 4, -4, -4, 11,
    -4, -4, 4, -11, 11, 4, -4, 4, 11, -4, 4, 4, -11, -11, -4, -4, -4, -11, -4,
    -4, -4, -11, 11, -4, -4, 4, -11, -4, 4, -4, -11,
};

/*
 * Gradients for 4D. They approximate the directions to the
 * vertices of a disprismatotesseractihexadecachoron from the center,
 * skewed so that the tetrahedral and cubic facets can be inscribed inside
 * spheres of the same radius.
 */
static const signed char gradients4D[] = {
    3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, -3, 1, 1, 1, -1, 3, 1, 1,
    -1, 1, 3, 1, -1, 1, 1, 3, 3, -1, 1, 1, 1, -3, 1, 1, 1, -1, 3, 1, 1, -1, 1,
    3, -3, -1, 1, 1, -1, -3, 1, 1, -1, -1, 3, 1, -1, -1, 1, 3, 3, 1, -1, 1, 1,
    3, -1, 1, 1, 1, -3, 1, 1, 1, -1, 3, -3, 1, -1, 1, -1, 3, -1, 1, -1, 1, -3,
    1, -1, 1, -1, 3, 3, -1, -1, 1, 1, -3, -1, 1, 1, -1, -3, 1, 1, -1, -1, 3, -3,
    -1, -1, 1, -1, -3, -1, 1, -1, -1, -3, 1, -1, -1, -1, 3, 3, 1, 1, -1, 1, 3,
    1, -1, 1, 1, 3, -1, 1, 1, 1, -3, -3, 1, 1, -1, -1, 3, 1, -1, -1, 1, 3, -1,
    -1, 1, 1, -3, 3, -1, 1, -1, 1, -3, 1, -1, 1, -1, 3, -1, 1, -1, 1, -3, -3,
    -1, 1, -1, -1, -3, 1, -1, -1, -1, 3, -1, -1, -1, 1, -3, 3, 1, -1, -1, 1, 3,
    -1, -1, 1, 1, -3, -1, 1, 1, -1, -3, -3, 1, -1, -1, -1, 3, -1, -1, -1, 1, -3,
    -1, -1, 1, -1, -3, 3, -1, -1, -1, 1, -3, -1, -1, 1, -1, -3, -1, 1, -1, -1,
    -3, -3, -1, -1, -1, -1, -3, -1, -1, -1, -1, -3, -1, -1, -1, -1, -3,
};

static double extrapolate2(
    struct osn_context* ctx, int xsb, int ysb, double dx, double dy)
{
    int16_t* perm = ctx->perm;
    int index = perm[(perm[xsb & 0xFF] + ysb) & 0xFF] & 0x0E;
    return gradients2D[index] * dx + gradients2D[index + 1] * dy;
}

static inline int fastFloor(double x)
{
    int xi = (int) x;
    return x < xi ? xi - 1 : xi;
}

static int allocate_perm(struct osn_context* ctx, int nperm, int ngrad)
{
    if (ctx->perm)
        free(ctx->perm);
    if (ctx->permGradIndex3D)
        free(ctx->permGradIndex3D);
    ctx->perm = (int16_t*) malloc(sizeof(*ctx->perm) * nperm);
    if (!ctx->perm)
        return -ENOMEM;
    ctx->permGradIndex3D =
        (int16_t*) malloc(sizeof(*ctx->permGradIndex3D) * ngrad);
    if (!ctx->permGradIndex3D) {
        free(ctx->perm);
        return -ENOMEM;
    }
    return 0;
}

static int par_simplex_noise(int64_t seed, struct osn_context** ctx)
{
    int rc;
    int16_t source[256];
    int i;
    int16_t* perm;
    int16_t* permGradIndex3D;

    *ctx = (struct osn_context*) malloc(sizeof(**ctx));
    if (!(*ctx))
        return -ENOMEM;
    (*ctx)->perm = NULL;
    (*ctx)->permGradIndex3D = NULL;

    rc = allocate_perm(*ctx, 256, 256);
    if (rc) {
        free(*ctx);
        return rc;
    }

    perm = (*ctx)->perm;
    permGradIndex3D = (*ctx)->permGradIndex3D;

    for (i = 0; i < 256; i++)
        source[i] = (int16_t) i;
    seed = seed * 6364136223846793005LL + 1442695040888963407LL;
    seed = seed * 6364136223846793005LL + 1442695040888963407LL;
    seed = seed * 6364136223846793005LL + 1442695040888963407LL;
    for (i = 255; i >= 0; i--) {
        seed = seed * 6364136223846793005LL + 1442695040888963407LL;
        int r = (int) ((seed + 31) % (i + 1));
        if (r < 0)
            r += (i + 1);
        perm[i] = source[r];
        permGradIndex3D[i] =
            (short) ((perm[i] % (ARRAYSIZE(gradients3D) / 3)) * 3);
        source[r] = source[i];
    }
    return 0;
}

static void par_simplex_noise_free(struct osn_context* ctx)
{
    if (!ctx)
        return;
    if (ctx->perm) {
        free(ctx->perm);
        ctx->perm = NULL;
    }
    if (ctx->permGradIndex3D) {
        free(ctx->permGradIndex3D);
        ctx->permGradIndex3D = NULL;
    }
    free(ctx);
}

static double par_simplex_noise2(struct osn_context* ctx, double x, double y)
{
    // Place input coordinates onto grid.
    double stretchOffset = (x + y) * STRETCH_CONSTANT_2D;
    double xs = x + stretchOffset;
    double ys = y + stretchOffset;

    // Floor to get grid coordinates of rhombus (stretched square) super-cell
    // origin.
    int xsb = fastFloor(xs);
    int ysb = fastFloor(ys);

    // Skew out to get actual coordinates of rhombus origin. We'll need these
    // later.
    double squishOffset = (xsb + ysb) * SQUISH_CONSTANT_2D;
    double xb = xsb + squishOffset;
    double yb = ysb + squishOffset;

    // Compute grid coordinates relative to rhombus origin.
    double xins = xs - xsb;
    double yins = ys - ysb;

    // Sum those together to get a value that determines which region we're in.
    double inSum = xins + yins;

    // Positions relative to origin point.
    double dx0 = x - xb;
    double dy0 = y - yb;

    // We'll be defining these inside the next block and using them afterwards.
    double dx_ext, dy_ext;
    int xsv_ext, ysv_ext;

    double value = 0;

    // Contribution (1,0)
    double dx1 = dx0 - 1 - SQUISH_CONSTANT_2D;
    double dy1 = dy0 - 0 - SQUISH_CONSTANT_2D;
    double attn1 = 2 - dx1 * dx1 - dy1 * dy1;
    if (attn1 > 0) {
        attn1 *= attn1;
        value += attn1 * attn1 * extrapolate2(ctx, xsb + 1, ysb + 0, dx1, dy1);
    }

    // Contribution (0,1)
    double dx2 = dx0 - 0 - SQUISH_CONSTANT_2D;
    double dy2 = dy0 - 1 - SQUISH_CONSTANT_2D;
    double attn2 = 2 - dx2 * dx2 - dy2 * dy2;
    if (attn2 > 0) {
        attn2 *= attn2;
        value += attn2 * attn2 * extrapolate2(ctx, xsb + 0, ysb + 1, dx2, dy2);
    }

    if (inSum <= 1) {  // We're inside the triangle (2-Simplex) at (0,0)
        double zins = 1 - inSum;
        if (zins > xins ||
            zins >
                yins) {  //(0,0) is one of the closest two triangular vertices
            if (xins > yins) {
                xsv_ext = xsb + 1;
                ysv_ext = ysb - 1;
                dx_ext = dx0 - 1;
                dy_ext = dy0 + 1;
            } else {
                xsv_ext = xsb - 1;
                ysv_ext = ysb + 1;
                dx_ext = dx0 + 1;
                dy_ext = dy0 - 1;
            }
        } else {  //(1,0) and (0,1) are the closest two vertices.
            xsv_ext = xsb + 1;
            ysv_ext = ysb + 1;
            dx_ext = dx0 - 1 - 2 * SQUISH_CONSTANT_2D;
            dy_ext = dy0 - 1 - 2 * SQUISH_CONSTANT_2D;
        }
    } else {  // We're inside the triangle (2-Simplex) at (1,1)
        double zins = 2 - inSum;
        if (zins < xins ||
            zins <
                yins) {  //(0,0) is one of the closest two triangular vertices
            if (xins > yins) {
                xsv_ext = xsb + 2;
                ysv_ext = ysb + 0;
                dx_ext = dx0 - 2 - 2 * SQUISH_CONSTANT_2D;
                dy_ext = dy0 + 0 - 2 * SQUISH_CONSTANT_2D;
            } else {
                xsv_ext = xsb + 0;
                ysv_ext = ysb + 2;
                dx_ext = dx0 + 0 - 2 * SQUISH_CONSTANT_2D;
                dy_ext = dy0 - 2 - 2 * SQUISH_CONSTANT_2D;
            }
        } else {  //(1,0) and (0,1) are the closest two vertices.
            dx_ext = dx0;
            dy_ext = dy0;
            xsv_ext = xsb;
            ysv_ext = ysb;
        }
        xsb += 1;
        ysb += 1;
        dx0 = dx0 - 1 - 2 * SQUISH_CONSTANT_2D;
        dy0 = dy0 - 1 - 2 * SQUISH_CONSTANT_2D;
    }

    // Contribution (0,0) or (1,1)
    double attn0 = 2 - dx0 * dx0 - dy0 * dy0;
    if (attn0 > 0) {
        attn0 *= attn0;
        value += attn0 * attn0 * extrapolate2(ctx, xsb, ysb, dx0, dy0);
    }

    // Extra Vertex
    double attn_ext = 2 - dx_ext * dx_ext - dy_ext * dy_ext;
    if (attn_ext > 0) {
        attn_ext *= attn_ext;
        value += attn_ext * attn_ext *
                 extrapolate2(ctx, xsv_ext, ysv_ext, dx_ext, dy_ext);
    }

    return value / NORM_CONSTANT_2D;
}

#undef PAR_MIN
#undef PAR_MAX
#undef PAR_CLAMP
#undef PAR_MALLOC
#undef PAR_CALLOC
#undef PAR_SWAP
#undef PAR_PI
#endif
