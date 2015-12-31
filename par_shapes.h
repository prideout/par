// SHAPES :: https://github.com/prideout/par
// Simple generation and transformation of indexed triangles.
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

#define PAR_shapes_SMOOTH_NORMALS (1 << 0)
#define PAR_shapes_FACET_NORMALS  (1 << 1)
#define PAR_shapes_TEXTURE_COORDS (1 << 2)

char const * const * par_shapes_list_parametric();
par_shapes_mesh const* par_shapes_create_parametric(char const*, int slices,
    int stacks, int flags);
void par_shapes_free(par_shapes_mesh const*);
void par_shapes_export(par_shapes_mesh const*, char const* objfile);

// Misc
par_shapes_mesh const* par_shapes_create_tree(int seed, int flags);
par_shapes_mesh const* par_shapes_create_rock(int seed, int flags);

// Transforms
void par_shapes_compute_normals(par_shapes_mesh*, int faceted);
void par_shapes_translate(par_shapes_mesh*, float x, float y, float z);
void par_shapes_scale(par_shapes_mesh*, float x, float y, float z);
void par_shapes_rotate(par_shapes_mesh*, float radians,
    float x, float y, float z);
void par_shapes_merge(par_shapes_mesh* dst, par_shapes_mesh const* src);

// http://prideout.net/blog/?p=44
par_shapes_mesh const* par_shapes_create_icosahedron();
par_shapes_mesh const* par_shapes_create_octohedron();
par_shapes_mesh const* par_shapes_create_cube();
par_shapes_mesh const* par_shapes_create_sphere(int nsubdivisions);
par_shapes_mesh const* par_shapes_subdivide(par_shapes_mesh const*);
// par_shapes_mesh const* par_shapes_create_tube_from_callback(...);
// par_shapes_mesh const* par_shapes_create_surf_from_callback(...);

// http://prideout.net/blog/?p=72
par_shapes_mesh const* par_shapes_create_lsystem(char const* program);

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

typedef void (*par_shapes_fn)(float* const, float*);

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

char const * const * par_shapes_list_parametric()
{
    return par_shapes_names;
}

par_shapes_mesh const* par_shapes_create_parametric(char const* name,
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

void par_shapes_free(par_shapes_mesh const* public_mesh)
{
    par_shapes_mesh* mesh = (par_shapes_mesh*) public_mesh;
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
    uint16_t const* indices = mesh->triangles;
    if (tcoords) {
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

#undef PAR_MIN
#undef PAR_MAX
#undef PAR_CLAMP
#undef PAR_ALLOC
#undef PAR_SWAP
#undef PAR_PI
#endif
