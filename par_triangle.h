// TRIANGLE :: https://github.com/prideout/par
// Consumes polygons, produces constrained Delaunay triangulations, etc.
//
// This implements "An improved incremental algorithm for constructing
// restricted Delaunay triangulations" by Marc Vigo Anglada.  We might
// eventually replace this with Shewchuk's randomized method, which is
// probably much faster.
//
// Here's an example that tessellates a pentagon with a triangular hole:
//
//    const uint16_t lengths[2] = {5, 3};
//    const float points[] = {
//        0, 0, 6, 0, 6, 4, 3, 7, 0, 4, // CCW pentagon
//        3, 3, 4, 1, 2, 1, // CW triangle
//    };
//    par_triangle_path* path = par_triangle_path_create(lengths, 2, points, 0);
//    par_triangle_mesh* mesh = par_triangle_mesh_create_cdt(path);
//    ...draw mesh here...
//    par_triangle_mesh_free(mesh);
//    par_triangle_path_free(path);
//
// TODO items
// --------------
// (0) CDT via Anglada's algorithm
//     https://www.dropbox.com/s/nkqsi1u7ey1bec7/Incremental_CDT.pdf
// (1) Refinement via Chew's Second Algorithm
//     https://en.wikipedia.org/wiki/Chew%27s_second_algorithm
// (2) Extrusion (and path reversal)
//     useful for msquares coastline (and a lighttransport/nanort demo)
// (3) Ear clipping
//     https://github.com/mapbox/earcut.hpp (uses boost)
//     https://github.com/prideout/polygon.js
// (4) Improve CDT generation via Shewchuk's randomized method
//     Fast Segment Insertion and Incremental Construction of
//     Constrained Delaunay Triangulations
//
// The MIT License
// Copyright (c) 2016 Philip Rideout

#ifndef PAR_TRIANGLE_H
#define PAR_TRIANGLE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

// Planar straight-line graph composed of one or more "loops" (polygons)
// whereby counterclockwise loops are solid and clockwise loops are holes.
// When serializing to SVG, all loops can be aggregated into a single <path>,
// provided they each terminate with "Z" and use the default fill rule.
typedef struct {
    float* points;       // list of XY vertex coordinates
    int npoints;         // number of 2-tuples in "points"
    float** loops;       // list of pointers to the start of each loop
    uint16_t* lengths;   // list of loop lengths
    int nloops;          // number of loops
} par_triangle_path;

typedef struct par_triangle_mesh_s {
    float const* points;       // Flat list of 2-tuples (X Y X Y...)
    int npoints;               // Number of points
    uint16_t const* triangles; // Flat list of 3-tuples (I J K I J K...)
    int ntriangles;            // Number of triangles
} par_triangle_mesh;

// Create a planar straight-line graph or resize an existing graph.  When
// creating a brand new graph, clients should pass null to "old".  The lengths
// argument is a client-owned array that specifies the desired number of
// two-tuples in each loop.  The provided points array can be null, in which
// case the client is expected to populate the coordinates after construction.
par_triangle_path* par_triangle_path_create(uint16_t const* lengths,
    int nloops, float const* points, par_triangle_path* old);

// Free all memory associated with a planar straight-line graph.
void par_triangle_path_free(par_triangle_path*);

// Triangulate the given polygon set using constrained Delaunay tessellation.
par_triangle_mesh* par_triangle_mesh_create_cdt(par_triangle_path const*);

// Free all memory associated with a 2D triangle mesh.
void par_triangle_mesh_free(par_triangle_mesh*);

// Scale, then translate, every point in the given path.
void par_triangle_path_transform(par_triangle_path* path, float sx, float sy,
    float tx, float ty);

// Find the triangle that contains the given point, otherwise return -1.
int par_triangle_mesh_find_triangle(par_triangle_mesh const* mesh, float x,
    float y);

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_TRIANGLE_IMPLEMENTATION
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <float.h>
#include <stdbool.h>

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

#ifndef pa_free
#define pa_free(a) ((a) ? PAR_FREE(pa___raw(a)), 0 : 0)
#define pa_push(a, v) (pa___maybegrow(a, 1), (a)[pa___n(a)++] = (v))
#define pa_pop(a) (pa___n(a)--)
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

par_triangle_path* par_triangle_path_create(uint16_t const* lengths,
    int nloops, float const* points, par_triangle_path* old)
{
    par_triangle_path* path = old;
    if (!path) {
        path = PAR_CALLOC(par_triangle_path, 1);
    }
    path->nloops = nloops;
    path->lengths = PAR_REALLOC(uint16_t, path->lengths, nloops);
    memcpy(path->lengths, lengths, sizeof(uint16_t) * nloops);
    path->npoints = 0;
    for (int i = 0; i < nloops; i++) {
        path->npoints += lengths[i];
    }
    path->points = PAR_REALLOC(float, path->points, path->npoints * 2);
    if (points) {
        memcpy(path->points, points, path->npoints * 2 * sizeof(float));
    }
    path->loops = PAR_REALLOC(float*, path->loops, nloops);
    float* pt = path->points;
    for (int i = 0; i < nloops; i++) {
        path->loops[i] = pt;
        pt += 2 * lengths[i];
    }
    return path;
}

void par_triangle_path_free(par_triangle_path* path)
{
    PAR_FREE(path->points);
    PAR_FREE(path->loops);
    PAR_FREE(path->lengths);
    PAR_FREE(path);
}

typedef struct par_triangle__vert_t {
    float x;
    float y;
    struct par_triangle__edge_t* outgoing;
} par_triangle__vert;

typedef struct par_triangle__face_t {
    struct par_triangle__edge_t* edge;
} par_triangle__face;

typedef struct par_triangle__edge_t {
    struct par_triangle__vert_t* end;
    struct par_triangle__face_t* face;
    struct par_triangle__edge_t* next;
    struct par_triangle__edge_t* pair;
} par_triangle__edge;

typedef struct {

    // Public data.
    float* points;
    int npoints;
    uint16_t* triangles;
    int ntriangles;

    // Private data; includes a half-edge data structure.
    par_triangle__vert* verts;
    par_triangle__face* faces;
    par_triangle__edge* edges;

} par_triangle__mesh;

// static float par_triangle__dot(float x1, float y1, float x2, float y2)
// {
//     return x1 * x2 + y1 * y2;
// }

static float par_triangle__cross(float x1, float y1, float x2, float y2)
{
    return x1 * y2 - y1 * x2;
}

// Assumes CCW ordering.
static bool par_triangle__contains(float px, float py, float x1, float y1,
    float x2, float y2, float x3, float y3)
{
    float vx1 = x2 - x1, vy1 = y2 - y1;
    float px1 = px - x1, py1 = py - y1;
    float vx2 = x3 - x2, vy2 = y3 - y2;
    float px2 = px - x2, py2 = py - y2;
    float vx3 = x1 - x3, vy3 = y1 - y3;
    float px3 = px - x3, py3 = py - y3;
    float c1 = par_triangle__cross(px1, py1, vx1, vy1);
    float c2 = par_triangle__cross(px2, py2, vx2, vy2);
    float c3 = par_triangle__cross(px3, py3, vx3, vy3);
    return c1 <= 0 && c2 <= 0 && c3 <= 0;
}

int par_triangle_mesh_find_triangle(par_triangle_mesh const* m, float x,
    float y)
{
    par_triangle__mesh const* mesh = (par_triangle__mesh const*) m;
    int nfaces = pa_count(mesh->faces);
    par_triangle__face const* face = mesh->faces;
    for (int f = 0; f < nfaces; f++, face++) {
        par_triangle__edge const* edge = face->edge;
        par_triangle__vert const* a = edge->end;
        par_triangle__vert const* b = edge->next->end;
        par_triangle__vert const* c = edge->next->next->end;
        if (par_triangle__contains(x, y, a->x, a->y, b->x, b->y, c->x, c->y)) {
            return f;
        }
    }
    return -1;
}

// Populate a brand new mesh with one triangle that wholly encompasses the unit
// square in [0, +1].
static par_triangle__mesh* par_triangle__mesh_create()
{
    par_triangle__mesh* result = PAR_CALLOC(par_triangle__mesh, 1);
    pa_add(result->verts, 3);
    pa_add(result->faces, 1);
    pa_add(result->edges, 3);
    const par_triangle__vert bigtri[] = {
        {-1.0, -1.0, 0},
        {+2.0, -1.0, 0},
        {+0.5, +3.0, 0}
    };
    memcpy(result->verts, bigtri, sizeof(bigtri));
    result->edges[0].end = result->verts + 1;
    result->edges[0].face = result->faces + 0;
    result->edges[0].next = result->edges + 1;
    result->edges[0].pair = 0;//result->edges + 3;
    result->edges[1].end = result->verts + 2;
    result->edges[1].face = result->faces + 0;
    result->edges[1].next = result->edges + 2;
    result->edges[1].pair = 0;//result->edges + 4;
    result->edges[2].end = result->verts + 0;
    result->edges[2].face = result->faces + 0;
    result->edges[2].next = result->edges + 0;
    result->edges[2].pair = 0;//result->edges + 5;
    result->faces[0].edge = result->edges + 0;
    result->verts[0].outgoing = result->edges + 0;
    result->verts[1].outgoing = result->edges + 1;
    result->verts[2].outgoing = result->edges + 2;
    // result->edges[3].end = result->verts + 0;
    // result->edges[3].face = 0;
    // result->edges[3].next = result->edges + 2;
    // result->edges[3].pair = result->edges + 0;
    // result->edges[4].end = result->verts + 1;
    // result->edges[4].face = 0;
    // result->edges[4].next = result->edges + 0;
    // result->edges[4].pair = result->edges + 1;
    // result->edges[5].end = result->verts + 2;
    // result->edges[5].face = 0;
    // result->edges[5].next = result->edges + 1;
    // result->edges[5].pair = result->edges + 2;
    return result;
}

static void par_triangle__mesh_validate(par_triangle__mesh* mesh)
{
    int nfaces = mesh->ntriangles = pa_count(mesh->faces);
    par_triangle__face const* face = mesh->faces;
    for (int f = 0; f < nfaces; f++, face++) {
        par_triangle__edge const* e = face->edge;
        par_triangle__vert const* a = e->end;
        par_triangle__vert const* b = e->next->end;
        par_triangle__vert const* c = e->next->next->end;
        float ab[2] = {b->x - a->x, b->y - a->y};
        float ac[2] = {c->x - a->x, c->y - a->y};
        assert(par_triangle__cross(ab[0], ab[1], ac[0], ac[1]) > 0);
        assert(!e->pair || e->pair->pair == e);
        assert(!e->next->pair || e->next->pair->pair == e->next);
        assert(!e->next->next->pair ||
            e->next->next->pair->pair == e->next->next);
        assert(e->face == face);
        assert(e->next->face == face);
        assert(e->next->next->face == face);
    }
    int nedges = pa_count(mesh->edges);
    int nborders = 0;
    par_triangle__edge const* edge = mesh->edges;
    for (int e = 0; e < nedges; e++, edge++) {
        if (!edge->pair) {
            nborders++;
            continue;
        }
        assert(edge->next->next->end == edge->pair->end);
    }
    assert(nborders == 3);
}

// Consume the half-edge data structure and produce data for the public fields.
static void par_triangle__mesh_finalize(par_triangle__mesh* mesh)
{
    // Produce "points" by consuming XY coordinates from the half-edge mesh.
    int nverts = mesh->npoints = pa_count(mesh->verts);
    float* point = mesh->points = PAR_MALLOC(float, 2 * nverts);
    par_triangle__vert const* vert = mesh->verts;
    for (int v = 0; v < nverts; v++, vert++, point += 2) {
        point[0] = vert->x;
        point[1] = vert->y;
    }

    // Produce "triangles" by consuming vertex pointers in the half-edge mesh.
    int nfaces = mesh->ntriangles = pa_count(mesh->faces);
    uint16_t* tri = mesh->triangles = PAR_MALLOC(uint16_t, 3 * nfaces);
    par_triangle__face const* face = mesh->faces;
    for (int f = 0; f < nfaces; f++, face++, tri += 3) {
        par_triangle__edge const* edge = face->edge;
        tri[0] = edge->end - mesh->verts;
        tri[1] = edge->next->end - mesh->verts;
        tri[2] = edge->next->next->end - mesh->verts;
    }
}

void par_triangle_path_transform(par_triangle_path* path, float sx, float sy,
    float tx, float ty)
{
    float* pt = path->points;
    for (int p = 0; p < path->npoints; p++, pt++) {
        pt[0] = pt[0] * sx + tx;
        pt[1] = pt[1] * sy + ty;
    }
}

static void par_triangle__mesh_transform(par_triangle__mesh* mesh,
    float sx, float sy, float tx, float ty)
{
    int nverts = pa_count(mesh->verts);
    par_triangle__vert* vert = mesh->verts;
    for (int v = 0; v < nverts; v++, vert++) {
        vert->x = vert->x * sx + tx;
        vert->y = vert->y * sy + ty;
    }
}

static void par_triangle__mesh_grow(par_triangle__mesh* mesh, int nverts,
    int nedges, int nfaces)
{
    // Reallocate verts and repair all pointers to verts.
    par_triangle__vert* verts = mesh->verts;
    pa_add(mesh->verts, nverts);
    if (verts != mesh->verts) {
        par_triangle__edge* edge = mesh->edges;
        for (int i = 0; i < pa_count(mesh->edges); i++, edge++) {
            edge->end = mesh->verts + (edge->end - verts);
        }
    }

    // Reallocate edges and repair all pointers to edges.
    par_triangle__edge* edges = mesh->edges;
    pa_add(mesh->edges, nedges);
    if (edges != mesh->edges) {
        par_triangle__edge* edge = mesh->edges;
        for (int i = 0; i < pa_count(mesh->edges) - nedges; i++, edge++) {
            edge->next = mesh->edges + (edge->next - edges);
            if (edge->pair) {
                edge->pair = mesh->edges + (edge->pair - edges);
            }
        }
        par_triangle__face* face = mesh->faces;
        for (int i = 0; i < pa_count(mesh->faces); i++, face++) {
            face->edge = mesh->edges + (face->edge - edges);
        }
        par_triangle__vert* vert = mesh->verts;
        for (int i = 0; i < pa_count(mesh->verts) - nverts; i++, vert++) {
            vert->outgoing = mesh->edges + (vert->outgoing - edges);
        }
    }

    // Reallocate faces and repair all pointers to faces.
    par_triangle__face* faces = mesh->faces;
    pa_add(mesh->faces, nfaces);
    if (faces != mesh->faces) {
        par_triangle__edge* edge = mesh->edges;
        for (int i = 0; i < pa_count(mesh->edges) - nedges; i++, edge++) {
            if (edge->face) {
                edge->face = mesh->faces + (edge->face - faces);
            }
        }
    }
}

// Change all edge pointers that were pointing to "from".
static void par_triangle__mesh_remap(par_triangle__edge* from,
    par_triangle__edge* to)
{
    if (from->next->next->end->outgoing == from) {
        from->next->next->end->outgoing = to;
    }
    if (from->pair) {
        from->pair->pair = to;
    }
    if (from->face && from->face->edge == from) {
        from->face->edge = to;
    }
    if (to) {
        *to = *from;
        from->next->next->next = to;
    }
}

// Remove a face and its three interior half-edges.
static void par_triangle__mesh_remove(par_triangle__mesh* mesh, int iface)
{
    int nedges = pa_count(mesh->edges);
    int nfaces = pa_count(mesh->faces);

    // Stash the edges that we're about to kill.
    par_triangle__edge* edgea0 = mesh->faces[iface].edge;
    par_triangle__edge* edgeb0 = edgea0->next;
    par_triangle__edge* edgec0 = edgeb0->next;

    // Stash the edges that we're about to move.
    par_triangle__edge* edgea1 = mesh->edges + nedges - 3;
    par_triangle__edge* edgeb1 = mesh->edges + nedges - 2;
    par_triangle__edge* edgec1 = mesh->edges + nedges - 1;

    // Move the last face into the slot currently occupied by the dead face.
    par_triangle__face* face0 = mesh->faces + iface;
    par_triangle__face* face1 = mesh->faces + nfaces - 1;
    face1->edge->face = face0;
    face1->edge->next->face = face0;
    face1->edge->next->next->face = face0;
    face0->edge = face1->edge;
    pa___n(mesh->faces) -= 1;

    // Remap all edge pointers appropriately.
    par_triangle__mesh_remap(edgea0, 0);
    par_triangle__mesh_remap(edgeb0, 0);
    par_triangle__mesh_remap(edgec0, 0);
    par_triangle__mesh_remap(edgea1, edgea0);
    par_triangle__mesh_remap(edgeb1, edgeb0);
    par_triangle__mesh_remap(edgec1, edgec0);
    pa___n(mesh->edges) -= 3;
}

static void par_triangle__mesh_subdivide(par_triangle__mesh* mesh, int face,
    float const* pt)
{
    // Stash the three vertices for the face that we're subdividing.
    par_triangle__edge* e = mesh->faces[face].edge;
    int av = e->end - mesh->verts;
    int bv = e->next->end - mesh->verts;
    int cv = e->next->next->end - mesh->verts;

    // Remove the face and its three half-edges.
    par_triangle__mesh_remove(mesh, face);

    // Add space for 1 new vertex, 9 new edges, and 3 new triangles.
    par_triangle__mesh_grow(mesh, 1, 9, 3);
    int nverts = pa_count(mesh->verts);
    int nedges = pa_count(mesh->edges);
    int nfaces = pa_count(mesh->faces);
    par_triangle__vert* vert = mesh->verts + nverts - 1;
    par_triangle__edge* edges = mesh->edges + nedges - 9;
    par_triangle__face* faces = mesh->faces + nfaces - 3;
    par_triangle__vert* v0 = mesh->verts + av;
    par_triangle__vert* v1 = mesh->verts + bv;
    par_triangle__vert* v2 = mesh->verts + cv;

    // New vertex.
    vert->x = pt[0];
    vert->y = pt[1];
    vert->outgoing = edges + 2;

    // New Face 0
    faces[0].edge = edges + 0;
    edges[0].face = faces + 0;
    edges[0].next = edges + 1;
    edges[0].end = v0;
    edges[1].face = faces + 0;
    edges[1].next = edges + 2;
    edges[1].end = vert;
    edges[2].face = faces + 0;
    edges[2].next = edges + 0;
    edges[2].end = v2;
    edges += 3;

    // New Face 1
    faces[1].edge = edges + 0;
    edges[0].face = faces + 1;
    edges[0].next = edges + 1;
    edges[0].end = v1;
    edges[1].face = faces + 1;
    edges[1].next = edges + 2;
    edges[1].end = vert;
    edges[2].face = faces + 1;
    edges[2].next = edges + 0;
    edges[2].end = v0;
    edges += 3;

    // New Face 2
    faces[2].edge = edges + 0;
    edges[0].face = faces + 2;
    edges[0].next = edges + 1;
    edges[0].end = v2;
    edges[1].face = faces + 2;
    edges[1].next = edges + 2;
    edges[1].end = vert;
    edges[2].face = faces + 2;
    edges[2].next = edges + 0;
    edges[2].end = v1;
    edges -= 6;

    // Populate the pair pointers.
    edges[0].pair = 0;
    edges[1].pair = edges + 5;
    edges[2].pair = edges + 7;
    edges[3].pair = 0;
    edges[4].pair = edges + 8;
    edges[5].pair = edges + 1;
    edges[6].pair = 0;
    edges[7].pair = edges + 2;
    edges[8].pair = edges + 4;
    par_triangle__edge* edge = mesh->edges;
    for (int i = 0; i < pa_count(mesh->edges) - 9; i++, edge++) {
        if (edge->end == v0 && edge->next->end == v2) {
            edge->next->pair = edges + 0;
            edges[0].pair = edge->next;
            continue;
        }
        if (edge->end == v1 && edge->next->end == v0) {
            edge->next->pair = edges + 3;
            edges[3].pair = edge->next;
            continue;
        }
        if (edge->end == v2 && edge->next->end == v1) {
            edge->next->pair = edges + 6;
            edges[6].pair = edge->next;
            continue;
        }
    }
}

// This is an implementation of Anglada's AddPointCDT function.
static void par_triangle__mesh_addpt(par_triangle__mesh* mesh, float const* pt)
{
    float x = pt[0];
    float y = pt[1];
    par_triangle_mesh* public_mesh = (par_triangle_mesh*) mesh;
    int tri = par_triangle_mesh_find_triangle(public_mesh, x, y);
    par_triangle__mesh_subdivide(mesh, tri, pt);
}

par_triangle_mesh* par_triangle_mesh_create_cdt(par_triangle_path const* path)
{
    par_triangle__mesh* mesh = par_triangle__mesh_create();
    par_triangle_mesh* result = (par_triangle_mesh*) mesh;
    float minx = FLT_MAX, maxx = -FLT_MAX;
    float miny = FLT_MAX, maxy = -FLT_MAX;
    float const* pt = path->points;
    for (int p = 0; p < path->npoints; p++, pt += 2) {
        minx = PAR_MIN(pt[0], minx);
        miny = PAR_MIN(pt[1], miny);
        maxx = PAR_MAX(pt[0], maxx);
        maxy = PAR_MAX(pt[1], maxy);
    }
    float width = maxx - minx;
    float height = maxy - miny;
    par_triangle__mesh_transform(mesh, width, height, minx, miny);
    pt = path->points;
    for (int p = 0; p < path->npoints; p++, pt += 2) {
        par_triangle__mesh_validate(mesh);
        par_triangle__mesh_addpt(mesh, pt);
    }
    par_triangle__mesh_finalize(mesh);
    return result;
}

void par_triangle_mesh_free(par_triangle_mesh* m)
{
    par_triangle__mesh* mesh = (par_triangle__mesh*) m;
    PAR_FREE(mesh->points);
    PAR_FREE(mesh->triangles);
    pa_free(mesh->verts);
    pa_free(mesh->faces);
    pa_free(mesh->edges);
    PAR_FREE(mesh);
}

#endif // PAR_TRIANGLE_IMPLEMENTATION
#endif // PAR_TRIANGLE_H
