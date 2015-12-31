// SHAPES :: https:github.com/prideout/par
// Mesh generator for parametric surfaces and other simple geometry.
//
//     http://github.prideout.net/c-shapes/
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

#include <stdint.h>

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

typedef struct {
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

// TBD, http://prideout.net/blog/?p=44
typedef void (*par_shapes_fn)(float* const, float*);
par_shapes_mesh* par_shapes_create_custom_parametric(par_shapes_fn, int slices,
    int stacks, int flags);
void par_shapes_compute_normals(par_shapes_mesh*);
par_shapes_mesh const* par_shapes_create_tree(int seed, int flags);
par_shapes_mesh const* par_shapes_create_rock(int seed, int flags);
par_shapes_mesh const* par_shapes_create_icosahedron();
par_shapes_mesh const* par_shapes_create_octohedron();
par_shapes_mesh const* par_shapes_create_cube(int inward); // for Cornell boxes
par_shapes_mesh const* par_shapes_create_sphere(int nsubdivisions);
par_shapes_mesh const* par_shapes_subdivide(par_shapes_mesh const*);

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_SHAPES_IMPLEMENTATION

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>

#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_ALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#define PAR_PI (3.14159265359)

static void par_shapes_private_sphere(float* const uv, float* xyz);
static void par_shapes_private_plane(float* const uv, float* xyz);
static void par_shapes_private_klein(float* const uv, float* xyz);
static void par_shapes_private_cylinder(float* const uv, float* xyz);
static void par_shapes_private_torus(float* const uv, float* xyz);

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

static void par_shapes_cross(float* result, float const* a, float const* b)
{
    result[0] = (a[1] * b[2]) - (a[2] * b[1]);
    result[1] = (a[2] * b[0]) - (a[0] * b[2]);
    result[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

static void par_shapes_normalize(float* v)
{
    float lsqr = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    if (lsqr > 0) {
        float scale = 1.0f / lsqr;
        v[0] *= scale;
        v[1] *= scale;
        v[2] *= scale;
    }
}

static void par_shapes_subtract(float* result, float const* a)
{
    result[0] -= a[0];
    result[1] -= a[1];
    result[2] -= a[2];
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
    par_shapes_mesh* mesh = calloc(sizeof(par_shapes_mesh), 1);

    // Generate verts.
    mesh->npoints = (slices + 1) * (stacks + 1);
    mesh->points = calloc(sizeof(float) * 3 * mesh->npoints, 1);
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
        mesh->normals = calloc(sizeof(float) * 3 * mesh->npoints, 1);
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
                par_shapes_subtract(du_xyz, xyz);
                par_shapes_subtract(dv_xyz, xyz);
                par_shapes_cross(normals, du_xyz, dv_xyz);
                par_shapes_normalize(normals);
                normals += 3;
            }
        }
    }

    // Generate texture coordinates.
    if (flags & PAR_SHAPES_TEXTURE_COORDS) {
        mesh->tcoords = calloc(sizeof(float) * 2 * mesh->npoints, 1);
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
    mesh->triangles = calloc(sizeof(uint16_t) * 3 * mesh->ntriangles, 1);
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
    float* points = malloc(vecsize * npoints);
    memcpy(points, dst->points, vecsize * dst->npoints);
    memcpy(points + 3 * dst->npoints, src->points, vecsize * src->npoints);
    free(dst->points);
    dst->points = points;
    dst->npoints = npoints;
    if (src->normals || dst->normals) {
        float* norms = calloc(vecsize * npoints, 1);
        if (dst->normals) {
            memcpy(norms, dst->normals, vecsize * dst->npoints);
        }
        if (src->normals) {
            memcpy(norms + 3 * dst->npoints, src->normals,
                vecsize * src->npoints);
        }
        free(dst->normals);
        dst->normals = norms;
    }
    if (src->tcoords || dst->tcoords) {
        int uvsize = sizeof(float) * 2;
        float* uvs = calloc(uvsize * npoints, 1);
        if (dst->tcoords) {
            memcpy(uvs, dst->tcoords, uvsize * dst->npoints);
        }
        if (src->tcoords) {
            memcpy(uvs + 2 * dst->npoints, src->tcoords, uvsize * src->npoints);
        }
        free(dst->tcoords);
        dst->tcoords = uvs;
    }
    int ntriangles = dst->ntriangles + src->ntriangles;
    int trisize = sizeof(uint16_t) * 3;
    uint16_t* triangles = malloc(trisize * ntriangles);
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
    par_shapes_mesh* mesh = calloc(sizeof(par_shapes_mesh), 1);
    mesh->npoints = slices + 1;
    mesh->points = malloc(sizeof(float) * 3 * mesh->npoints);
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
    par_shapes_normalize(nnormal);
    if (flags & PAR_SHAPES_SMOOTH_NORMALS) {
        mesh->normals = malloc(sizeof(float) * 3 * mesh->npoints);
        float* norms = mesh->normals;
        for (int i = 0; i < mesh->npoints; i++) {
            *norms++ = nnormal[0];
            *norms++ = nnormal[1];
            *norms++ = nnormal[2];
        }
    }
    mesh->ntriangles = slices;
    mesh->triangles = malloc(sizeof(uint16_t) * 3 * mesh->ntriangles);
    uint16_t* triangles = mesh->triangles;
    for (int i = 0; i < slices; i++) {
        *triangles++ = 0;
        *triangles++ = 1 + i;
        *triangles++ = 1 + (i + 1) % slices;
    }
    float k[3] = {0, 0, -1};
    float axis[3];
    par_shapes_cross(axis, nnormal, k);
    par_shapes_normalize(axis);
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

#undef PAR_MIN
#undef PAR_MAX
#undef PAR_CLAMP
#undef PAR_ALLOC
#undef PAR_SWAP
#undef PAR_PI
#endif
