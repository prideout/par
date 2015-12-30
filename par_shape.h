// SHAPES :: https://github.com/prideout/par
// Simple generation and manipulation of indexed triangles.
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
} par_shape_mesh;

#define PAR_SHAPE_SMOOTH_NORMALS (1 << 0)
#define PAR_SHAPE_FACET_NORMALS (1 << 1)
#define PAR_SHAPE_TEXTURE_COORDS (1 << 2)

// http://prideout.net/blog/?p=44
void par_shape_free(par_shape_mesh const*);
const char** par_shape_list_parametric();
par_shape_mesh const* par_shape_create_parametric(char const*, int slices,
    int stacks, int flags);

// Misc
par_shape_mesh const* par_shape_create_tree(int seed, int flags);
par_shape_mesh const* par_shape_create_rock(int seed, int flags);

// Modifiers
void par_shape_compute_normals(par_shape_mesh*, int faceted);
void par_shape_translate(par_shape_mesh*, float x, float y, float z);
void par_shape_scale(par_shape_mesh*, float x, float y, float z);
void par_shape_rotate(par_shape_mesh*, float radians,
    float x, float y, float z);

// http://prideout.net/blog/?p=44
par_shape_mesh const* par_shape_create_icosahedron();
par_shape_mesh const* par_shape_create_octohedron();
par_shape_mesh const* par_shape_create_cube();
par_shape_mesh const* par_shape_create_sphere(int nsubdivisions);
par_shape_mesh const* par_shape_subdivide(par_shape_mesh const*);
// par_shape_mesh const* par_shape_create_tube_from_callback(...);
// par_shape_mesh const* par_shape_create_surf_from_callback(...);

// http://prideout.net/blog/?p=72
par_shape_mesh const* par_shape_create_lsystem(char const* program);

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_SHAPES_IMPLEMENTATION

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <string.h>

#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_ALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }

static const char* par_shape_shapes[] = {
    "sphere",
    "disk",
    "klein",
    "cylinder",
    "torus",
};


#undef PAR_MIN
#undef PAR_MAX
#undef PAR_CLAMP
#undef PAR_ALLOC
#undef PAR_SWAP
#endif
