// BUBBLES :: https://github.com/prideout/par
// Simple C library for packing circles into hierarchical (or flat) diagrams.
//
// Based on "Visualization of Large Hierarchical Data by Circle Packing" by
// Wang et al, with an extension for cylinderical space.
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

// Tiny POD structure returned by all packing functions.  Private data is
// attached after the public fields, so clients should call the provided
// free function rather than freeing the memory manually.
typedef struct {
    double* xyr; // array of 3-tuples (x y radius) in same order as input data
    int count;   // number of 3-tuples
} par_bubbles_t;

void par_bubbles_free_result(par_bubbles_t*);

// Flat Packing ----------------------------------------------------------------

// Entry point for simple non-hierarchical packing.  Takes a list of radii.
par_bubbles_t* par_bubbles_pack(double const* radiuses, int nradiuses);

// Hierarchical Packing --------------------------------------------------------

// All these functions consume a list of nodes, where each node is
// represented by an integer.  The integer is an index to the node's parent.
// Unlike the flat packing API, clients do not have control over radius.
// However, they do have control over bounding area.

// Pack the circles into a circle that is centered at the origin.
par_bubbles_t par_bubbles_hpack_circle(int* nodes, int nnodes, double radius);

// Pack the circles into a rectangle centered at the origin.  Takes a pointer
// to two doubles: width and height.
par_bubbles_t par_bubbles_hpack_rect(int* nodes, int nnodes, double* dims);

// Pack the circles into an unwrapped cylinder surface.  The dimensions are a
// two-tuple of circumference (i.e., unwrapped width) and height.
par_bubbles_t par_bubbles_hpack_cylinder(int* nodes, int nnodes, double* dims);

// Queries ---------------------------------------------------------------------

// Find the node at the given position.  Children are "on top" of their parents.
// If the result is -1, there is no node at the given pick coordinate.
int par_bubbles_pick(par_bubbles_t*, double x, double y);

// Take a pointer to 6 floats and set them to min xyz, max xyz.
void par_bubbles_compute_aabb(par_bubbles_t const*, double* aabb);

// Dump out a SVG file for diagnostic purposes.
void par_bubbles_export(par_bubbles_t const* bubbles, char const* filename);

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_BUBBLES_IMPLEMENTATION

#include <math.h>
#include <stdio.h>

typedef struct {
    int prev;
    int next;
} par_bubbles__node;

typedef struct {
    double* xyr;
    int count;
    double const* client_radii;
    int const* client_parents;
    par_bubbles__node* chain;
    int npacked;
} par_bubbles__t;

#ifndef PARB_HELPERS
#define PARB_HELPERS 1
#define PAR_PI (3.14159265359)
#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_MALLOC(T, N) ((T*) malloc(N * sizeof(T)))
#define PAR_CALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))
#define PAR_REALLOC(T, BUF, N) ((T*) realloc(BUF, sizeof(T) * N))
#define PAR_FREE(BUF) free(BUF)
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#endif

static void par_bubbles__place(double* c, double const* a, double const* b)
{
    double db = a[2] + c[2], dx = b[0] - a[0], dy = b[1] - a[1];
    if (db && (dx || dy)) {
        double da = b[2] + c[2], dc = dx * dx + dy * dy;
        da *= da;
        db *= db;
        double x = 0.5 + (db - da) / (2 * dc);
        double db1 = db - dc;
        double y0 = PAR_MAX(0, 2 * da * (db + dc) - db1 * db1 - da * da);
        double y = sqrt(y0) / (2 * dc);
        c[0] = a[0] + x * dx + y * dy;
        c[1] = a[1] + x * dy - y * dx;
    } else {
        c[0] = a[0] + db;
        c[1] = a[1];
    }
}

// static double par_bubbles__dist2(double const* a, double const* b)
// {
//     return (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]);
// }

static double par_bubbles__len2(double const* a)
{
    return a[0] * a[0] + a[1] * a[1];
}

static void par_bubbles__init(par_bubbles__t* bubbles)
{
    double* xyr = bubbles->xyr;
    double const* radii = bubbles->client_radii;
    par_bubbles__node* chain = bubbles->chain;

    *xyr++ = -*radii;
    *xyr++ = 0;
    *xyr++ = *radii++;
    if (bubbles->count == ++bubbles->npacked) {
        return;
    }

    *xyr++ = *radii;
    *xyr++ = 0;
    *xyr++ = *radii++;
    if (bubbles->count == ++bubbles->npacked) {
        return;
    }

    xyr[2] = *radii;
    par_bubbles__place(xyr, xyr - 6, xyr - 3);
    if (bubbles->count == ++bubbles->npacked) {
        return;
    }

    chain[0].prev = 2;
    chain[0].next = 1;
    chain[1].prev = 0;
    chain[1].next = 2;
    chain[2].prev = 1;
    chain[2].next = 0;
}

static void par_bubbles__pack(par_bubbles__t* bubbles)
{
    double* xyr = bubbles->xyr;
    double const* radii = bubbles->client_radii;
    par_bubbles__node* chain = bubbles->chain;
    int chain_length = 3;
    double dist, mindist;

    // Find the circle closest to the origin, known as "Cm" in the paper.
    int cm = 0;
    mindist = par_bubbles__len2(xyr + 0);
    dist = par_bubbles__len2(xyr + 3);
    if (dist > mindist) {
        cm = 1;
    }
    dist = par_bubbles__len2(xyr + 6);
    if (dist > mindist) {
        cm = 2;
    }

    // In the paper, "Cn" is always the node that follows Cm.
    int cn = chain[cm].next;

    for (; bubbles->npacked < bubbles->count; bubbles->npacked++) {
        // TODO
    }

    bubbles->count = 3; // TODO
}

par_bubbles_t* par_bubbles_pack(double const* radiuses, int nradiuses)
{
    par_bubbles__t* bubbles = PAR_CALLOC(par_bubbles__t, 1);
    if (nradiuses > 0) {
        bubbles->client_radii = radiuses;
        bubbles->count = nradiuses;
        bubbles->chain = PAR_MALLOC(par_bubbles__node, nradiuses);
        bubbles->xyr = PAR_MALLOC(double, 3 * nradiuses);
        par_bubbles__init(bubbles);
        par_bubbles__pack(bubbles);
    }
    return (par_bubbles_t*) bubbles;
}

void par_bubbles_free_result(par_bubbles_t* pubbub)
{
    par_bubbles__t* bubbles = (par_bubbles__t*) pubbub;
    PAR_FREE(bubbles->chain);
    PAR_FREE(bubbles->xyr);
    PAR_FREE(bubbles);
}

void par_bubbles_compute_aabb(par_bubbles_t const* bubbles, double* aabb)
{
    if (bubbles->count == 0) {
        return;
    }
    double* xyr = bubbles->xyr;
    aabb[0] = aabb[2] = xyr[0];
    aabb[1] = aabb[3] = xyr[1];
    for (int i = 0; i < bubbles->count; i++, xyr += 3) {
        aabb[0] = PAR_MIN(xyr[0] - xyr[2], aabb[0]);
        aabb[1] = PAR_MIN(xyr[1] - xyr[2], aabb[1]);
        aabb[2] = PAR_MAX(xyr[0] + xyr[2], aabb[2]);
        aabb[3] = PAR_MAX(xyr[1] + xyr[2], aabb[3]);
    }
}

void par_bubbles_export(par_bubbles_t const* bubbles, char const* filename)
{
    double aabb[4];
    par_bubbles_compute_aabb(bubbles, aabb);
    double maxextent = PAR_MAX(aabb[2] - aabb[0], aabb[3] - aabb[1]);
    double padding = 0.05 * maxextent;
    double spacing = 1.0 / maxextent;
    FILE* svgfile = fopen(filename, "wt");
    fprintf(svgfile,
        "<svg viewBox='%f %f %f %f' width='500px' height='500px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g stroke-width='%f' stroke='white' fill-opacity='0.2' fill='white'>\n"
        "<rect fill-opacity='1' stroke='#475473' fill='#475473' x='%f' y='%f' "
        "width='100%%' height='100%%'/>\n",
        aabb[0] - padding, aabb[1] - padding,
        aabb[2] - aabb[0] + 2 * padding, aabb[3] - aabb[1] + 2 * padding,
        spacing,
        aabb[0] - padding, aabb[1] - padding);
    double const* xyr = bubbles->xyr;
    for (int i = 0; i < bubbles->count; i++, xyr += 3) {
        fprintf(svgfile, "<circle cx='%f' cy='%f' r='%f'/>\n",
            xyr[0], xyr[1], xyr[2] - spacing);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
}

#endif
