// BUBBLES :: https://github.com/prideout/par
// Simple C library for packing circles into hierarchical (or flat) diagrams.
//
// Based on "Visualization of Large Hierarchical Data by Circle Packing" by
// Wang et al (2006), with an extension for cylinderical space.
//
// Also uses Emo Welzl's "Smallest enclosing disks" algorithm (1991).
//
// The API is divided into three sections:
//
//   - Enclosing.  Compute the smallest bounding circle for points or circles.
//   - Packing.    Pack circles together, or into other circles.
//   - Queries.    Given a touch points, pick a circle from a hierarchy, etc.
//
// In addition to the comment block above each function declaration, the API
// has informal documentation here:
//
//     http://github.prideout.net/bubbles/
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

#ifndef PAR_BUBBLES_H
#define PAR_BUBBLES_H

// Enclosing -------------------------------------------------------------------

// Reads an array of (x,y) coordinates, writes a single 3-tuple (x,y,radius).
void par_bubbles_enclose_points(double const* xy, int npts, double* result);

// Reads an array of 3-tuples (x,y,radius), writes a 3-tuple (x,y,radius).
void par_bubbles_enclose_disks(double const* xyr, int ndisks, double* result);

// Low-level function used internally to create a circle tangent to 3 points.
void par_bubbles_touch_three_points(double const* xy, double* result);

// Packing ---------------------------------------------------------------------

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
// The root node is its own parent, and it must be the first node in the list.
// Unlike the flat packing API, clients do not have control over radius.
// However, they do have control over bounding area.

// Pack the circles into a circle that is centered at the origin.
par_bubbles_t* par_bubbles_hpack_circle(int* nodes, int nnodes, double radius);

// Pack the circles into a rectangle centered at the origin.  Takes a pointer
// to two doubles: width and height.
par_bubbles_t* par_bubbles_hpack_rect(int* nodes, int nnodes, double* dims);

// Pack the circles into an unwrapped cylinder surface.  The dimensions are a
// two-tuple of circumference (i.e., unwrapped width) and height.
par_bubbles_t* par_bubbles_hpack_cylinder(int* nodes, int nnodes, double* dims);

// Queries ---------------------------------------------------------------------

// Find the node at the given position.  Children are "on top" of their parents.
// If the result is -1, there is no node at the given pick coordinate.
int par_bubbles_pick(par_bubbles_t*, double x, double y);

// Get bounding box; take a pointer to 4 floats and set them to min xy, max xy.
void par_bubbles_compute_aabb(par_bubbles_t const*, double* aabb);

// Dump out a SVG file for diagnostic purposes.
void par_bubbles_export(par_bubbles_t const* bubbles, char const* filename);

#ifndef PAR_HELPERS
#define PAR_HELPERS 1
#define PAR_PI (3.14159265359)
#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_MALLOC(T, N) ((T*) malloc(N * sizeof(T)))
#define PAR_CALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))
#define PAR_REALLOC(T, BUF, N) ((T*) realloc(BUF, sizeof(T) * N))
#define PAR_FREE(BUF) free(BUF)
#define PAR_SWAP(T, A, B) { T tmp = B; B = A; A = tmp; }
#define PAR_SQR(a) (a * a)
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#endif // PAR_BUBBLES_H

#ifdef PAR_BUBBLES_IMPLEMENTATION

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int prev;
    int next;
} par_bubbles__node;

typedef struct {
    double* xyr;              // results array
    int count;                // client-provided count
    double const* radiuses;   // client-provided radius list
    par_bubbles__node* chain; // counterclockwise enveloping chain
    int const* graph_parents; // client-provided parent indices
    int* graph_children;      // flat list of children indices
    int* graph_heads;         // list of "pointers" to first child
    int* graph_tails;         // list of "pointers" to one-past-last child
    int npacked;
    int maxwidth;
} par_bubbles__t;

// Assigns an xy to "c" such that it becomes tangent to "a" and "b".
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

static double par_bubbles__len2(double const* a)
{
    return a[0] * a[0] + a[1] * a[1];
}

static void par_bubbles__initgraph(par_bubbles__t* bubbles)
{
    int const* parents = bubbles->graph_parents;
    int* nchildren = PAR_CALLOC(int, bubbles->count);
    for (int i = 0; i < bubbles->count; i++) {
        nchildren[parents[i]]++;
    }
    int c = 0;
    bubbles->graph_heads = PAR_CALLOC(int, bubbles->count * 2);
    bubbles->graph_tails = bubbles->graph_heads + bubbles->count;
    for (int i = 0; i < bubbles->count; i++) {
        bubbles->maxwidth = PAR_MAX(bubbles->maxwidth, nchildren[i]);
        bubbles->graph_heads[i] = bubbles->graph_tails[i] = c;
        c += nchildren[i];
    }
    bubbles->graph_heads[0] = bubbles->graph_tails[0] = 1;
    bubbles->graph_children = PAR_MALLOC(int, c);
    for (int i = 1; i < bubbles->count; i++) {
        int parent = parents[i];
        bubbles->graph_children[bubbles->graph_tails[parent]++] = i;
    }
    PAR_FREE(nchildren);
}

static void par_bubbles__initflat(par_bubbles__t* bubbles)
{
    double* xyr = bubbles->xyr;
    double const* radii = bubbles->radiuses;
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

// March forward or backward along the enveloping chain, starting with the
// node at "cn" and testing for collision against the node at "ci".
static int par_bubbles__collide(par_bubbles__t* bubbles, int ci, int cn,
    int* cj, int direction)
{
    double const* ci_xyr = bubbles->xyr + ci * 3;
    par_bubbles__node* chain = bubbles->chain;
    int nsteps = 1;
    if (direction > 0) {
        for (int i = chain[cn].next; i != cn; i = chain[i].next, ++nsteps) {
            double const* i_xyr = bubbles->xyr + i * 3;
            double dx = i_xyr[0] - ci_xyr[0];
            double dy = i_xyr[1] - ci_xyr[1];
            double dr = i_xyr[2] + ci_xyr[2];
            if (0.999 * dr * dr > dx * dx + dy * dy) {
                *cj = i;
                return nsteps;
            }
        }
        return 0;
    }
    for (int i = chain[cn].prev; i != cn; i = chain[i].prev, ++nsteps) {
        double const* i_xyr = bubbles->xyr + i * 3;
        double dx = i_xyr[0] - ci_xyr[0];
        double dy = i_xyr[1] - ci_xyr[1];
        double dr = i_xyr[2] + ci_xyr[2];
        if (0.999 * dr * dr > dx * dx + dy * dy) {
            *cj = i;
            return nsteps;
        }
    }
    return 0;
}

static void par_bubbles__packflat(par_bubbles__t* bubbles)
{
    double const* radii = bubbles->radiuses;
    double* xyr = bubbles->xyr;
    par_bubbles__node* chain = bubbles->chain;

    // Find the circle closest to the origin, known as "Cm" in the paper.
    int cm = 0;
    double mindist = par_bubbles__len2(xyr + 0);
    double dist = par_bubbles__len2(xyr + 3);
    if (dist > mindist) {
        cm = 1;
    }
    dist = par_bubbles__len2(xyr + 6);
    if (dist > mindist) {
        cm = 2;
    }

    // In the paper, "Cn" is always the node that follows "Cm".
    int ci, cn = chain[cm].next;

    for (ci = bubbles->npacked; ci < bubbles->count; ) {
        double* ci_xyr = xyr + ci * 3;
        ci_xyr[2] = radii[ci];
        double* cm_xyr = xyr + cm * 3;
        double* cn_xyr = xyr + cn * 3;
        par_bubbles__place(ci_xyr, cn_xyr, cm_xyr);

        // Check for a collision.  In the paper, "Cj" is the intersecting node.
        int cj_f;
        int nfsteps = par_bubbles__collide(bubbles, ci, cn, &cj_f, +1);
        if (!nfsteps) {
            chain[cm].next = ci;
            chain[ci].prev = cm;
            chain[ci].next = cn;
            chain[cn].prev = ci;
            cm = ci++;
            continue;
        }

        // Search backwards for a collision, in case it is closer.
        int cj_b;
        int nbsteps = par_bubbles__collide(bubbles, ci, cm, &cj_b, -1);

        // Intersection occurred after Cn.
        if (nfsteps <= nbsteps) {
            cn = cj_f;
            chain[cm].next = cn;
            chain[cn].prev = cm;
            continue;
        }

        // Intersection occurred before Cm.
        cm = cj_b;
        chain[cm].next = cn;
        chain[cn].prev = cm;
    }

    bubbles->npacked = bubbles->count;
}

void par_bubbles_touch_three_points(double const* xy, double* xyr)
{
    // Many thanks to Stephen Schmitts:
    // http://www.abecedarical.com/zenosamples/zs_circle3pts.html
    double p1x = xy[0], p1y = xy[1];
    double p2x = xy[2], p2y = xy[3];
    double p3x = xy[4], p3y = xy[5];
    double a = p2x - p1x, b = p2y - p1y;
    double c = p3x - p1x, d = p3y - p1y;
    double e = a * (p2x + p1x) * 0.5 + b * (p2y + p1y) * 0.5;
    double f = c * (p3x + p1x) * 0.5 + d * (p3y + p1y) * 0.5;
    double det = a*d - b*c;
    double cx = xyr[0] = (d*e - b*f) / det;
    double cy = xyr[1] = (-c*e + a*f) / det;
    xyr[2] = sqrt((p1x - cx)*(p1x - cx) + (p1y - cy)*(p1y - cy));
}

par_bubbles_t* par_bubbles_pack(double const* radiuses, int nradiuses)
{
    par_bubbles__t* bubbles = PAR_CALLOC(par_bubbles__t, 1);
    if (nradiuses > 0) {
        bubbles->radiuses = radiuses;
        bubbles->count = nradiuses;
        bubbles->chain = PAR_MALLOC(par_bubbles__node, nradiuses);
        bubbles->xyr = PAR_MALLOC(double, 3 * nradiuses);
        par_bubbles__initflat(bubbles);
        par_bubbles__packflat(bubbles);
    }
    return (par_bubbles_t*) bubbles;
}

// TODO: use a stack instead of recursion
void par_bubbles__compute_radius(par_bubbles__t* bubbles,
    par_bubbles__t* worker, int parent)
{
    int head = bubbles->graph_heads[parent];
    int tail = bubbles->graph_tails[parent];
    int nchildren = tail - head;
    int pr = parent * 3 + 2;
    bubbles->xyr[pr] = 1;
    if (nchildren == 0) {
        return;
    }
    for (int children_index = head; children_index != tail; children_index++) {
        int child = bubbles->graph_children[children_index];
        par_bubbles__compute_radius(bubbles, worker, child);
        bubbles->xyr[pr] += bubbles->xyr[child * 3 + 2];
    }
    bubbles->xyr[pr] = sqrtf(bubbles->xyr[pr]);
}

void par_bubbles__hpack(par_bubbles__t* bubbles, par_bubbles__t* worker,
    int parent)
{
    int head = bubbles->graph_heads[parent];
    int tail = bubbles->graph_tails[parent];
    int nchildren = tail - head;
    if (nchildren == 0) {
        return;
    }

    // Cast away const because we're using the worker as a cache to avoid
    // a kazillion malloc / free calls.
    double* radiuses = (double*) worker->radiuses;

    // We perform flat layout twice: once without padding (to determine scale)
    // and then again with scaled padding.
    double cx, cy, cr;
    double px = bubbles->xyr[parent * 3 + 0];
    double py = bubbles->xyr[parent * 3 + 1];
    double pr = bubbles->xyr[parent * 3 + 2];
    const double PAR_HPACK_PADDING1 = 0.1;
    const double PAR_HPACK_PADDING2 = 0.05;
    double scaled_padding = 0.0;
    while (1) {
        worker->npacked = 0;
        worker->count = nchildren;
        int c = 0;
        for (int cindex = head; cindex != tail; cindex++) {
            int child = bubbles->graph_children[cindex];
            radiuses[c++] = bubbles->xyr[child * 3 + 2] + scaled_padding;
        }
        par_bubbles__initflat(worker);
        par_bubbles__packflat(worker);
        double aabb[6];
        par_bubbles_compute_aabb((par_bubbles_t const*) worker, aabb);
        cx = 0.5 * (aabb[0] + aabb[2]);
        cy = 0.5 * (aabb[1] + aabb[3]);
        cr = 0;
        for (int c = 0; c < nchildren; c++) {
            double x = worker->xyr[c * 3 + 0] - cx;
            double y = worker->xyr[c * 3 + 1] - cy;
            double r = worker->xyr[c * 3 + 2];
            cr = PAR_MAX(cr, r + sqrtf(x * x + y * y));
        }
        if (scaled_padding || !PAR_HPACK_PADDING1) {
            break;
        } else {
            scaled_padding = PAR_HPACK_PADDING1 / cr;
        }
    }
    scaled_padding *= cr;
    cr += PAR_HPACK_PADDING2 * cr;
    if (nchildren == 1) {
        cr *= 2;
    }

    // Transform the children to fit nicely into their parent.
    int c = 0;
    pr /= cr;
    for (int cindex = head; cindex != tail; cindex++) {
        int child = bubbles->graph_children[cindex];
        bubbles->xyr[child * 3 + 0] = px + pr * (worker->xyr[c * 3 + 0] - cx);
        bubbles->xyr[child * 3 + 1] = py + pr * (worker->xyr[c * 3 + 1] - cy);
        bubbles->xyr[child * 3 + 2] -= scaled_padding;
        bubbles->xyr[child * 3 + 2] *= pr;
        c++;
    }

    // Recursion.  TODO: It might be better to use our own stack here.
    for (int cindex = head; cindex != tail; cindex++) {
        par_bubbles__hpack(bubbles, worker, bubbles->graph_children[cindex]);
    }
}

par_bubbles_t* par_bubbles_hpack_circle(int* nodes, int nnodes, double radius)
{
    par_bubbles__t* bubbles = PAR_CALLOC(par_bubbles__t, 1);
    if (nnodes > 0) {
        bubbles->graph_parents = nodes;
        bubbles->count = nnodes;
        bubbles->chain = PAR_MALLOC(par_bubbles__node, nnodes);
        bubbles->xyr = PAR_MALLOC(double, 3 * nnodes);
        par_bubbles__initgraph(bubbles);
        par_bubbles__t* worker = PAR_CALLOC(par_bubbles__t, 1);
        worker->radiuses = PAR_MALLOC(double, bubbles->maxwidth);
        worker->chain = PAR_MALLOC(par_bubbles__node, bubbles->maxwidth);
        worker->xyr = PAR_MALLOC(double, 3 * bubbles->maxwidth);
        par_bubbles__compute_radius(bubbles, worker, 0);
        bubbles->xyr[0] = 0;
        bubbles->xyr[1] = 0;
        bubbles->xyr[2] = radius;
        par_bubbles__hpack(bubbles, worker, 0);
        par_bubbles_free_result((par_bubbles_t*) worker);
    }
    return (par_bubbles_t*) bubbles;
}

void par_bubbles_free_result(par_bubbles_t* pubbub)
{
    par_bubbles__t* bubbles = (par_bubbles__t*) pubbub;
    PAR_FREE(bubbles->graph_children);
    PAR_FREE(bubbles->graph_heads);
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
    FILE* svgfile = fopen(filename, "wt");
    fprintf(svgfile,
        "<svg viewBox='%f %f %f %f' width='700px' height='700px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g stroke-width='0.5' stroke-opacity='0.5' stroke='black' "
        "fill-opacity='0.2' fill='#2A8BB6' "
        "transform='scale(1 -1) transform(0 %f)'>\n"
        "<rect fill-opacity='0.1' stroke='none' fill='#2A8BB6' x='%f' y='%f' "
        "width='100%%' height='100%%'/>\n",
        aabb[0] - padding, aabb[1] - padding,
        aabb[2] - aabb[0] + 2 * padding, aabb[3] - aabb[1] + 2 * padding,
        aabb[1] - padding,
        aabb[0] - padding, aabb[1] - padding);
    double const* xyr = bubbles->xyr;
    for (int i = 0; i < bubbles->count; i++, xyr += 3) {
        fprintf(svgfile, "<circle stroke-width='%f' cx='%f' cy='%f' r='%f'/>\n",
            xyr[2] * 0.01, xyr[0], xyr[1], xyr[2]);
        fprintf(svgfile, "<text text-anchor='middle' stroke='none' "
            "x='%f' y='%f' font-size='%f'>%d</text>\n",
            xyr[0], xyr[1] + xyr[2] * 0.125, xyr[2] * 0.5, i);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
}

#endif // PAR_BUBBLES_IMPLEMENTATION
