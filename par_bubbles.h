// BUBBLES :: https://github.com/prideout/par
// Simple C library for packing circles into hierarchical (or flat) diagrams.
//
// Based on "Visualization of Large Hierarchical Data by Circle Packing" by
// Wang et al (2006).
//
// Also implements Emo Welzl's "Smallest enclosing disks" algorithm (1991).
//
// The API is divided into three sections:
//
//   - Enclosing.  Compute the smallest bounding circle for points or circles.
//   - Packing.    Pack circles together, or into other circles.
//   - Queries.    Given a touch point, pick a circle from a hierarchy, etc.
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

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

// This can be any signed integer type.
#ifndef PAR_BUBBLES_INT
#define PAR_BUBBLES_INT int32_t
#endif

// This must be "float" or "double" or "long double".
#ifndef PAR_BUBBLES_FLT
#define PAR_BUBBLES_FLT double
#endif

// Enclosing / Touching --------------------------------------------------------

// Read an array of (x,y) coordinates, write a single 3-tuple (x,y,radius).
void par_bubbles_enclose_points(PAR_BUBBLES_FLT const* xy, PAR_BUBBLES_INT npts,
    PAR_BUBBLES_FLT* result);

// Read an array of 3-tuples (x,y,radius), write a 3-tuple (x,y,radius).
// Internally, this approximates each disk with an enclosing octagon.
void par_bubbles_enclose_disks(PAR_BUBBLES_FLT const* xyr,
    PAR_BUBBLES_INT ndisks, PAR_BUBBLES_FLT* result);

// Find the circle (x,y,radius) that is tangent to 3 points (x,y).
void par_bubbles_touch_three_points(PAR_BUBBLES_FLT const* xy,
    PAR_BUBBLES_FLT* result);

// Find a position for disk "c" that makes it tangent to "a" and "b".
// Note that the ordering of a and b can affect where c will land.
// All three arguments are pointers to three-tuples (x,y,radius).
void par_bubbles_touch_two_disks(PAR_BUBBLES_FLT* c, PAR_BUBBLES_FLT const* a,
    PAR_BUBBLES_FLT const* b);

// Packing ---------------------------------------------------------------------

// Tiny POD structure returned by all packing functions.  Private data is
// attached after the public fields, so clients should call the provided
// free function rather than freeing the memory manually.
typedef struct {
    PAR_BUBBLES_FLT* xyr;  // array of 3-tuples (x y radius) in input order
    PAR_BUBBLES_INT count; // number of 3-tuples in "xyr"
    PAR_BUBBLES_INT* ids;  // populated by par_bubbles_cull
} par_bubbles_t;

void par_bubbles_free_result(par_bubbles_t*);

// Entry point for unbounded non-hierarchical packing.  Takes a list of radii.
par_bubbles_t* par_bubbles_pack(PAR_BUBBLES_FLT const* radiuses,
    PAR_BUBBLES_INT nradiuses);

// Consume a hierarchy defined by a list of integers.  Each integer is an index
// to its parent. The root node is its own parent, and it must be the first node
// in the list. Clients do not have control over individual radiuses, only the
// radius of the outermost enclosing disk.
par_bubbles_t* par_bubbles_hpack_circle(PAR_BUBBLES_INT* nodes,
    PAR_BUBBLES_INT nnodes, PAR_BUBBLES_FLT radius);

// Queries ---------------------------------------------------------------------

// Find the node at the given position.  Children are on top of their parents.
// If the result is -1, there is no node at the given pick coordinate.
PAR_BUBBLES_INT par_bubbles_pick(par_bubbles_t const*, PAR_BUBBLES_FLT x,
    PAR_BUBBLES_FLT y);

// Get bounding box; take a pointer to 4 floats and set them to min xy, max xy.
void par_bubbles_compute_aabb(par_bubbles_t const*, PAR_BUBBLES_FLT* aabb);

// Check if the given circle (3-tuple) intersects the given aabb (4-tuple).
bool par_bubbles_check_aabb(PAR_BUBBLES_FLT const* disk,
    PAR_BUBBLES_FLT const* aabb);

// Clip the bubble diagram to the given AABB (4-tuple of left,bottom,right,top)
// and return the result.  Circles smaller than the given world-space
// "minradius" are removed.  Optionally, an existing diagram (dst) can be passed
// in to receive the culled dataset, which reduces the number of memory allocs
// when calling this function frequently.  Pass null to "dst" to create a new
// culled diagram.
par_bubbles_t* par_bubbles_cull(par_bubbles_t const* src,
    PAR_BUBBLES_FLT const* aabb, PAR_BUBBLES_FLT minradius, par_bubbles_t* dst);

// Dump out a SVG file for diagnostic purposes.
void par_bubbles_export(par_bubbles_t const* bubbles, char const* filename);
void par_bubbles_export_local(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT idx, char const* filename);

// Returns a pointer to a list of children nodes.
void par_bubbles_get_children(par_bubbles_t const* bubbles, PAR_BUBBLES_INT idx,
    PAR_BUBBLES_INT** pchildren, PAR_BUBBLES_INT* nchildren);

// Finds the height of the tree and returns one of its deepest leaves.
void par_bubbles_get_maxdepth(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT* maxdepth, PAR_BUBBLES_INT* leaf);

// Finds the height of the tree at a certain node.
PAR_BUBBLES_INT par_bubbles_get_depth(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT node);

// Returns a 4-tuple (min xy, max xy) for the given node.
void par_bubbles_compute_aabb_for_node(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT node, PAR_BUBBLES_FLT* aabb);

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
#define PAR_REALLOC(T, BUF, N) ((T*) realloc(BUF, sizeof(T) * N))
#define PAR_FREE(BUF) free(BUF)
#endif

#ifdef __cplusplus
}
#endif
#endif // PAR_BUBBLES_H

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_BUBBLES_IMPLEMENTATION
#define PARINT PAR_BUBBLES_INT
#define PARFLT PAR_BUBBLES_FLT

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>

typedef struct {
    PARINT prev;
    PARINT next;
} par_bubbles__node;

typedef struct {
    PARFLT* xyr;              // results array
    PARINT count;                // client-provided count
    PARINT* ids;                 // populated by par_bubbles_cull
    PARFLT const* radiuses;   // client-provided radius list
    par_bubbles__node* chain; // counterclockwise enveloping chain
    PARINT const* graph_parents; // client-provided parent indices
    PARINT* graph_children;      // flat list of children indices
    PARINT* graph_heads;         // list of "pointers" to first child
    PARINT* graph_tails;         // list of "pointers" to one-past-last child
    PARINT npacked;
    PARINT maxwidth;
    PARINT capacity;
} par_bubbles__t;

static PARFLT par_bubbles__len2(PARFLT const* a)
{
    return a[0] * a[0] + a[1] * a[1];
}

static void par_bubbles__initgraph(par_bubbles__t* bubbles)
{
    PARINT const* parents = bubbles->graph_parents;
    PARINT* nchildren = PAR_CALLOC(PARINT, bubbles->count);
    for (PARINT i = 0; i < bubbles->count; i++) {
        nchildren[parents[i]]++;
    }
    PARINT c = 0;
    bubbles->graph_heads = PAR_CALLOC(PARINT, bubbles->count * 2);
    bubbles->graph_tails = bubbles->graph_heads + bubbles->count;
    for (PARINT i = 0; i < bubbles->count; i++) {
        bubbles->maxwidth = PAR_MAX(bubbles->maxwidth, nchildren[i]);
        bubbles->graph_heads[i] = bubbles->graph_tails[i] = c;
        c += nchildren[i];
    }
    bubbles->graph_heads[0] = bubbles->graph_tails[0] = 1;
    bubbles->graph_children = PAR_MALLOC(PARINT, c);
    for (PARINT i = 1; i < bubbles->count; i++) {
        PARINT parent = parents[i];
        bubbles->graph_children[bubbles->graph_tails[parent]++] = i;
    }
    PAR_FREE(nchildren);
}

static void par_bubbles__initflat(par_bubbles__t* bubbles)
{
    PARFLT* xyr = bubbles->xyr;
    PARFLT const* radii = bubbles->radiuses;
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
    par_bubbles_touch_two_disks(xyr, xyr - 6, xyr - 3);
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
static PARINT par_bubbles__collide(par_bubbles__t* bubbles, PARINT ci,
    PARINT cn, PARINT* cj, PARINT direction)
{
    PARFLT const* ci_xyr = bubbles->xyr + ci * 3;
    par_bubbles__node* chain = bubbles->chain;
    PARINT nsteps = 1;
    if (direction > 0) {
        for (PARINT i = chain[cn].next; i != cn; i = chain[i].next, ++nsteps) {
            PARFLT const* i_xyr = bubbles->xyr + i * 3;
            PARFLT dx = i_xyr[0] - ci_xyr[0];
            PARFLT dy = i_xyr[1] - ci_xyr[1];
            PARFLT dr = i_xyr[2] + ci_xyr[2];
            if (0.999 * dr * dr > dx * dx + dy * dy) {
                *cj = i;
                return nsteps;
            }
        }
        return 0;
    }
    for (PARINT i = chain[cn].prev; i != cn; i = chain[i].prev, ++nsteps) {
        PARFLT const* i_xyr = bubbles->xyr + i * 3;
        PARFLT dx = i_xyr[0] - ci_xyr[0];
        PARFLT dy = i_xyr[1] - ci_xyr[1];
        PARFLT dr = i_xyr[2] + ci_xyr[2];
        if (0.999 * dr * dr > dx * dx + dy * dy) {
            *cj = i;
            return nsteps;
        }
    }
    return 0;
}

static void par_bubbles__packflat(par_bubbles__t* bubbles)
{
    PARFLT const* radii = bubbles->radiuses;
    PARFLT* xyr = bubbles->xyr;
    par_bubbles__node* chain = bubbles->chain;

    // Find the circle closest to the origin, known as "Cm" in the paper.
    PARINT cm = 0;
    PARFLT mindist = par_bubbles__len2(xyr + 0);
    PARFLT dist = par_bubbles__len2(xyr + 3);
    if (dist > mindist) {
        cm = 1;
    }
    dist = par_bubbles__len2(xyr + 6);
    if (dist > mindist) {
        cm = 2;
    }

    // In the paper, "Cn" is always the node that follows "Cm".
    PARINT ci, cn = chain[cm].next;

    for (ci = bubbles->npacked; ci < bubbles->count; ) {
        PARFLT* ci_xyr = xyr + ci * 3;
        ci_xyr[2] = radii[ci];
        PARFLT* cm_xyr = xyr + cm * 3;
        PARFLT* cn_xyr = xyr + cn * 3;
        par_bubbles_touch_two_disks(ci_xyr, cn_xyr, cm_xyr);

        // Check for a collision.  In the paper, "Cj" is the intersecting node.
        PARINT cj_f;
        PARINT nfsteps = par_bubbles__collide(bubbles, ci, cn, &cj_f, +1);
        if (!nfsteps) {
            chain[cm].next = ci;
            chain[ci].prev = cm;
            chain[ci].next = cn;
            chain[cn].prev = ci;
            cm = ci++;
            continue;
        }

        // Search backwards for a collision, in case it is closer.
        PARINT cj_b;
        PARINT nbsteps = par_bubbles__collide(bubbles, ci, cm, &cj_b, -1);

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

static void par__disk_from_two(PARFLT const* xy1, PARFLT const* xy2,
    PARFLT* result)
{
    PARFLT dx = xy1[0] - xy2[0];
    PARFLT dy = xy1[1] - xy2[1];
    result[0] = 0.5 * (xy1[0] + xy2[0]);
    result[1] = 0.5 * (xy1[1] + xy2[1]);
    result[2] = sqrt(dx * dx + dy * dy) / 2.0;
}

static PARINT par__disk_contains(PARFLT const* xyr, PARFLT const* xy)
{
    PARFLT dx = xyr[0] - xy[0];
    PARFLT dy = xyr[1] - xy[1];
    return dx * dx + dy * dy <= PAR_SQR(xyr[2]);
}

static void par__easydisk(PARFLT* disk, PARFLT const* edgepts, PARINT nedgepts)
{
    if (nedgepts == 0) {
        disk[0] = 0;
        disk[1] = 0;
        disk[2] = 0;
        return;
    }
    if (nedgepts == 1) {
        disk[0] = edgepts[0];
        disk[1] = edgepts[1];
        disk[2] = 0;
        return;
    }
    par__disk_from_two(edgepts, edgepts + 2, disk);
    if (nedgepts == 2 || par__disk_contains(disk, edgepts + 4)) {
        return;
    }
    par__disk_from_two(edgepts, edgepts + 4, disk);
    if (par__disk_contains(disk, edgepts + 2)) {
        return;
    }
    par__disk_from_two(edgepts + 2, edgepts + 4, disk);
    if (par__disk_contains(disk, edgepts)) {
        return;
    }
    par_bubbles_touch_three_points(edgepts, disk);
}

static void par__minidisk(PARFLT* disk, PARFLT const* pts, PARINT npts,
    PARFLT const* edgepts, PARINT nedgepts)
{
    if (npts == 0 || nedgepts == 3) {
        par__easydisk(disk, edgepts, nedgepts);
        return;
    }
    PARFLT const* pt = pts + (--npts) * 2;
    par__minidisk(disk, pts, npts, edgepts, nedgepts);
    if (!par__disk_contains(disk, pt)) {
        PARFLT edgepts1[6];
        for (PARINT i = 0; i < nedgepts * 2; i += 2) {
            edgepts1[i] = edgepts[i];
            edgepts1[i + 1] = edgepts[i + 1];
        }
        edgepts1[2 * nedgepts] = pt[0];
        edgepts1[2 * nedgepts + 1] = pt[1];
        par__minidisk(disk, pts, npts, edgepts1, ++nedgepts);
    }
}

static void par_bubbles__copy_disk(par_bubbles__t const* src,
    par_bubbles__t* dst, PARINT parent)
{
    PARINT i = dst->count++;
    if (dst->capacity < dst->count) {
        dst->capacity = PAR_MAX(16, dst->capacity) * 2;
        dst->xyr = PAR_REALLOC(PARFLT, dst->xyr, 3 * dst->capacity);
        dst->ids = PAR_REALLOC(PARINT, dst->ids, dst->capacity);
    }
    PARFLT const* xyr = src->xyr + parent * 3;
    dst->xyr[i * 3] = xyr[0];
    dst->xyr[i * 3 + 1] = xyr[1];
    dst->xyr[i * 3 + 2] = xyr[2];
    dst->ids[i] = parent;
}

void par_bubbles_enclose_points(PARFLT const* xy, PARINT npts, PARFLT* result)
{
    if (npts == 0) {
        return;
    }
    par__minidisk(result, xy, npts, 0, 0);
}

void par_bubbles_enclose_disks(PARFLT const* xyr, PARINT ndisks, PARFLT* result)
{
    PARINT ngon = 8;
    PARINT npts = ndisks * ngon;
    PARFLT* pts = PAR_MALLOC(PARFLT, npts * 2);
    PARFLT* ppts = pts;
    float dtheta = PAR_PI * 2.0 / ngon;

    for (PARINT i = 0; i < ndisks; i++) {
        PARFLT cx = xyr[i * 3];
        PARFLT cy = xyr[i * 3 + 1];
        PARFLT cr = xyr[i * 3 + 2];
        PARFLT a = 2.0 * cr / (1.0 + sqrt(2));
        PARFLT r = 0.5 * sqrt(2) * a * sqrt(2 + sqrt(2));
        float theta = 0;
        for (PARINT j = 0; j < ngon; j++, theta += dtheta) {
            *ppts++ = cx + r * cos(theta);
            *ppts++ = cy + r * sin(theta);
        }
    }
    par_bubbles_enclose_points(pts, npts, result);
    PAR_FREE(pts);
}

void par_bubbles_touch_three_points(PARFLT const* xy, PARFLT* xyr)
{
    // Many thanks to Stephen Schmitts:
    // http://www.abecedarical.com/zenosamples/zs_circle3pts.html
    PARFLT p1x = xy[0], p1y = xy[1];
    PARFLT p2x = xy[2], p2y = xy[3];
    PARFLT p3x = xy[4], p3y = xy[5];
    PARFLT a = p2x - p1x, b = p2y - p1y;
    PARFLT c = p3x - p1x, d = p3y - p1y;
    PARFLT e = a * (p2x + p1x) * 0.5 + b * (p2y + p1y) * 0.5;
    PARFLT f = c * (p3x + p1x) * 0.5 + d * (p3y + p1y) * 0.5;
    PARFLT det = a*d - b*c;
    PARFLT cx = xyr[0] = (d*e - b*f) / det;
    PARFLT cy = xyr[1] = (-c*e + a*f) / det;
    xyr[2] = sqrt((p1x - cx)*(p1x - cx) + (p1y - cy)*(p1y - cy));
}

void par_bubbles_touch_two_disks(PARFLT* c, PARFLT const* a, PARFLT const* b)
{
    PARFLT db = a[2] + c[2], dx = b[0] - a[0], dy = b[1] - a[1];
    if (db && (dx || dy)) {
        PARFLT da = b[2] + c[2], dc = dx * dx + dy * dy;
        da *= da;
        db *= db;
        PARFLT x = 0.5 + (db - da) / (2 * dc);
        PARFLT db1 = db - dc;
        PARFLT y0 = PAR_MAX(0, 2 * da * (db + dc) - db1 * db1 - da * da);
        PARFLT y = sqrt(y0) / (2 * dc);
        c[0] = a[0] + x * dx + y * dy;
        c[1] = a[1] + x * dy - y * dx;
    } else {
        c[0] = a[0] + db;
        c[1] = a[1];
    }
}

void par_bubbles_free_result(par_bubbles_t* pubbub)
{
    par_bubbles__t* bubbles = (par_bubbles__t*) pubbub;
    PAR_FREE(bubbles->graph_children);
    PAR_FREE(bubbles->graph_heads);
    PAR_FREE(bubbles->chain);
    PAR_FREE(bubbles->xyr);
    PAR_FREE(bubbles->ids);
    PAR_FREE(bubbles);
}

par_bubbles_t* par_bubbles_pack(PARFLT const* radiuses, PARINT nradiuses)
{
    par_bubbles__t* bubbles = PAR_CALLOC(par_bubbles__t, 1);
    if (nradiuses > 0) {
        bubbles->radiuses = radiuses;
        bubbles->count = nradiuses;
        bubbles->chain = PAR_MALLOC(par_bubbles__node, nradiuses);
        bubbles->xyr = PAR_MALLOC(PARFLT, 3 * nradiuses);
        par_bubbles__initflat(bubbles);
        par_bubbles__packflat(bubbles);
    }
    return (par_bubbles_t*) bubbles;
}

// TODO: use a stack instead of recursion
void par_bubbles__compute_radius(par_bubbles__t* bubbles,
    par_bubbles__t* worker, PARINT parent)
{
    PARINT head = bubbles->graph_heads[parent];
    PARINT tail = bubbles->graph_tails[parent];
    PARINT nchildren = tail - head;
    PARINT pr = parent * 3 + 2;
    bubbles->xyr[pr] = 1;
    if (nchildren == 0) {
        return;
    }
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = bubbles->graph_children[cindex];
        par_bubbles__compute_radius(bubbles, worker, child);
        bubbles->xyr[pr] += bubbles->xyr[child * 3 + 2];
    }
    bubbles->xyr[pr] = sqrtf(bubbles->xyr[pr]);
}

void par_bubbles__hpack(par_bubbles__t* bubbles, par_bubbles__t* worker,
    PARINT parent)
{
    PARINT head = bubbles->graph_heads[parent];
    PARINT tail = bubbles->graph_tails[parent];
    PARINT nchildren = tail - head;
    if (nchildren == 0) {
        return;
    }

    // Cast away const because we're using the worker as a cache to avoid
    // a kazillion malloc / free calls.
    PARFLT* radiuses = (PARFLT*) worker->radiuses;

    // We perform flat layout twice: once without padding (to determine scale)
    // and then again with scaled padding.
    PARFLT enclosure[3];
    PARFLT px = bubbles->xyr[parent * 3 + 0];
    PARFLT py = bubbles->xyr[parent * 3 + 1];
    PARFLT pr = bubbles->xyr[parent * 3 + 2];
    const PARFLT PAR_HPACK_PADDING1 = 0.1;
    const PARFLT PAR_HPACK_PADDING2 = 0.05;
    PARFLT scaled_padding = 0.0;
    while (1) {
        worker->npacked = 0;
        worker->count = nchildren;
        PARINT c = 0;
        for (PARINT cindex = head; cindex != tail; cindex++) {
            PARINT child = bubbles->graph_children[cindex];
            radiuses[c++] = bubbles->xyr[child * 3 + 2] + scaled_padding;
        }
        par_bubbles__initflat(worker);
        par_bubbles__packflat(worker);

        // Using Welzl's algorithm instead of a simple AABB enclosure is
        // slightly slower and doesn't yield much aesthetic improvement.

        #if PAR_BUBBLES_HPACK_WELZL
        par_bubbles_enclose_disks(worker->xyr, nchildren, enclosure);
        #else
        PARFLT aabb[6];
        par_bubbles_compute_aabb((par_bubbles_t const*) worker, aabb);
        enclosure[0] = 0.5 * (aabb[0] + aabb[2]);
        enclosure[1] = 0.5 * (aabb[1] + aabb[3]);
        enclosure[2] = 0;
        for (PARINT c = 0; c < nchildren; c++) {
            PARFLT x = worker->xyr[c * 3 + 0] - enclosure[0];
            PARFLT y = worker->xyr[c * 3 + 1] - enclosure[1];
            PARFLT r = worker->xyr[c * 3 + 2];
            enclosure[2] = PAR_MAX(enclosure[2], r + sqrtf(x * x + y * y));
        }
        #endif

        if (scaled_padding || !PAR_HPACK_PADDING1) {
            break;
        } else {
            scaled_padding = PAR_HPACK_PADDING1 / enclosure[2];
        }
    }
    PARFLT cx = enclosure[0], cy = enclosure[1], cr = enclosure[2];
    scaled_padding *= cr;
    cr += PAR_HPACK_PADDING2 * cr;
    if (nchildren == 1) {
        cr *= 2;
    }

    // Transform the children to fit nicely into their parent.
    PARINT c = 0;
    pr /= cr;
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = bubbles->graph_children[cindex];
        bubbles->xyr[child * 3 + 0] = px + pr * (worker->xyr[c * 3 + 0] - cx);
        bubbles->xyr[child * 3 + 1] = py + pr * (worker->xyr[c * 3 + 1] - cy);
        bubbles->xyr[child * 3 + 2] -= scaled_padding;
        bubbles->xyr[child * 3 + 2] *= pr;
        c++;
    }

    // Recursion.  TODO: It might be better to use our own stack here.
    for (PARINT cindex = head; cindex != tail; cindex++) {
        par_bubbles__hpack(bubbles, worker, bubbles->graph_children[cindex]);
    }
}

par_bubbles_t* par_bubbles_hpack_circle(PARINT* nodes, PARINT nnodes,
    PARFLT radius)
{
    par_bubbles__t* bubbles = PAR_CALLOC(par_bubbles__t, 1);
    if (nnodes > 0) {
        bubbles->graph_parents = nodes;
        bubbles->count = nnodes;
        bubbles->chain = PAR_MALLOC(par_bubbles__node, nnodes);
        bubbles->xyr = PAR_MALLOC(PARFLT, 3 * nnodes);
        par_bubbles__initgraph(bubbles);
        par_bubbles__t* worker = PAR_CALLOC(par_bubbles__t, 1);
        worker->radiuses = PAR_MALLOC(PARFLT, bubbles->maxwidth);
        worker->chain = PAR_MALLOC(par_bubbles__node, bubbles->maxwidth);
        worker->xyr = PAR_MALLOC(PARFLT, 3 * bubbles->maxwidth);
        par_bubbles__compute_radius(bubbles, worker, 0);
        bubbles->xyr[0] = 0;
        bubbles->xyr[1] = 0;
        bubbles->xyr[2] = radius;
        par_bubbles__hpack(bubbles, worker, 0);
        par_bubbles_free_result((par_bubbles_t*) worker);
    }
    return (par_bubbles_t*) bubbles;
}

// TODO: use a stack instead of recursion
static PARINT par_bubbles__pick(par_bubbles__t const* bubbles, PARINT parent,
    PARFLT x, PARFLT y)
{
    PARFLT const* xyr = bubbles->xyr + parent * 3;
    PARFLT d2 = PAR_SQR(x - xyr[0]) + PAR_SQR(y - xyr[1]);
    if (d2 > PAR_SQR(xyr[2])) {
        return -1;
    }
    PARINT head = bubbles->graph_heads[parent];
    PARINT tail = bubbles->graph_tails[parent];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = bubbles->graph_children[cindex];
        PARINT result = par_bubbles__pick(bubbles, child, x, y);
        if (result > -1) {
            return result;
        }
    }
    return parent;
}

PARINT par_bubbles_pick(par_bubbles_t const* cbubbles, PARFLT x, PARFLT y)
{
    par_bubbles__t const* bubbles = (par_bubbles__t const*) cbubbles;
    if (bubbles->count == 0) {
        return -1;
    }
    return par_bubbles__pick(bubbles, 0, x, y);
}

void par_bubbles_compute_aabb(par_bubbles_t const* bubbles, PARFLT* aabb)
{
    if (bubbles->count == 0) {
        return;
    }
    PARFLT const* xyr = bubbles->xyr;
    aabb[0] = aabb[2] = xyr[0];
    aabb[1] = aabb[3] = xyr[1];
    for (PARINT i = 0; i < bubbles->count; i++, xyr += 3) {
        aabb[0] = PAR_MIN(xyr[0] - xyr[2], aabb[0]);
        aabb[1] = PAR_MIN(xyr[1] - xyr[2], aabb[1]);
        aabb[2] = PAR_MAX(xyr[0] + xyr[2], aabb[2]);
        aabb[3] = PAR_MAX(xyr[1] + xyr[2], aabb[3]);
    }
}

bool par_bubbles_check_aabb(PARFLT const* disk, PARFLT const* aabb)
{
    PARFLT cx = PAR_CLAMP(disk[0], aabb[0], aabb[2]);
    PARFLT cy = PAR_CLAMP(disk[1], aabb[1], aabb[3]);
    PARFLT dx = disk[0] - cx;
    PARFLT dy = disk[1] - cy;
    PARFLT d2 = dx * dx + dy * dy;
    return d2 < (disk[2] * disk[2]);
}

static void par_bubbles__cull(par_bubbles__t const* src, PARFLT const* aabb,
    PARFLT minradius, par_bubbles__t* dst, PARINT parent)
{
    PARFLT const* xyr = src->xyr + parent * 3;
    if (xyr[2] < minradius || !par_bubbles_check_aabb(xyr, aabb)) {
        return;
    }
    par_bubbles__copy_disk(src, dst, parent);
    PARINT head = src->graph_heads[parent];
    PARINT tail = src->graph_tails[parent];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        par_bubbles__cull(src, aabb, minradius, dst, child);
    }
}

par_bubbles_t* par_bubbles_cull(par_bubbles_t const* psrc,
    PARFLT const* aabb, PARFLT minradius, par_bubbles_t* pdst)
{
    par_bubbles__t const* src = (par_bubbles__t const*) psrc;
    par_bubbles__t* dst = (par_bubbles__t*) pdst;
    if (!dst) {
        dst = PAR_CALLOC(par_bubbles__t, 1);
        pdst = (par_bubbles_t*) dst;
    } else {
        dst->count = 0;
    }
    if (src->count == 0) {
        return pdst;
    }
    par_bubbles__cull(src, aabb, minradius, dst, 0);
    return pdst;
}

void par_bubbles_export(par_bubbles_t const* bubbles, char const* filename)
{
    PARFLT aabb[4];
    par_bubbles_compute_aabb(bubbles, aabb);
    PARFLT maxextent = PAR_MAX(aabb[2] - aabb[0], aabb[3] - aabb[1]);
    PARFLT padding = 0.05 * maxextent;
    FILE* svgfile = fopen(filename, "wt");
    fprintf(svgfile,
        "<svg viewBox='%f %f %f %f' width='700px' height='700px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g stroke-width='0.5' stroke-opacity='0.5' stroke='black' "
        "fill-opacity='0.2' fill='#2A8BB6'>\n"
        "<rect fill-opacity='0.1' stroke='none' fill='#2A8BB6' x='%f' y='%f' "
        "width='100%%' height='100%%'/>\n",
        aabb[0] - padding, aabb[1] - padding,
        aabb[2] - aabb[0] + 2 * padding, aabb[3] - aabb[1] + 2 * padding,
        aabb[0] - padding, aabb[1] - padding);
    PARFLT const* xyr = bubbles->xyr;
    for (PARINT i = 0; i < bubbles->count; i++, xyr += 3) {
        fprintf(svgfile, "<circle stroke-width='%f' cx='%f' cy='%f' r='%f'/>\n",
            xyr[2] * 0.01, xyr[0], xyr[1], xyr[2]);
        fprintf(svgfile, "<text text-anchor='middle' stroke='none' "
            "x='%f' y='%f' font-size='%f'>%d</text>\n",
            xyr[0], xyr[1] + xyr[2] * 0.125, xyr[2] * 0.5, (int) i);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
}

void par_bubbles_get_children(par_bubbles_t const* pbubbles, PARINT node,
    PARINT** pchildren, PARINT* nchildren)
{
    par_bubbles__t const* bubbles = (par_bubbles__t const*) pbubbles;
    *pchildren = bubbles->graph_children + bubbles->graph_heads[node];
    *nchildren = bubbles->graph_tails[node] - bubbles->graph_heads[node];
}

void par_bubbles__get_maxdepth(par_bubbles__t const* bubbles, PARINT* maxdepth,
    PARINT* leaf, PARINT parent, PARINT depth)
{
    if (depth > *maxdepth) {
        *leaf = parent;
        *maxdepth = depth;
    }
    PARINT* children;
    PARINT nchildren;
    par_bubbles_t const* pbubbles = (par_bubbles_t const*) bubbles;
    par_bubbles_get_children(pbubbles, parent, &children, &nchildren);
    for (PARINT c = 0; c < nchildren; c++) {
        par_bubbles__get_maxdepth(bubbles, maxdepth, leaf, children[c],
            depth + 1);
    }
}

void par_bubbles_get_maxdepth(par_bubbles_t const* pbubbles, PARINT* maxdepth,
    PARINT* leaf)
{
    par_bubbles__t const* bubbles = (par_bubbles__t const*) pbubbles;
    *maxdepth = -1;
    *leaf = -1;
    return par_bubbles__get_maxdepth(bubbles, maxdepth, leaf, 0, 0);
}

PARINT par_bubbles_get_depth(par_bubbles_t const* pbubbles, PARINT node)
{
    par_bubbles__t const* bubbles = (par_bubbles__t const*) pbubbles;
    PARINT const* parents = bubbles->graph_parents;
    PARINT depth = 0;
    while (node) {
        node = parents[node];
        depth++;
    }
    return depth;
}

void par_bubbles_compute_aabb_for_node(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT node, PAR_BUBBLES_FLT* aabb)
{
    PARFLT const* xyr = bubbles->xyr + 3 * node;
    aabb[0] = aabb[2] = xyr[0];
    aabb[1] = aabb[3] = xyr[1];
    aabb[0] = PAR_MIN(xyr[0] - xyr[2], aabb[0]);
    aabb[1] = PAR_MIN(xyr[1] - xyr[2], aabb[1]);
    aabb[2] = PAR_MAX(xyr[0] + xyr[2], aabb[2]);
    aabb[3] = PAR_MAX(xyr[1] + xyr[2], aabb[3]);
}

void par_bubbles_export_local(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT idx, char const* filename)
{
    PAR_BUBBLES_FLT aabb[4];
    par_bubbles_compute_aabb_for_node(bubbles, idx, aabb);
    PARFLT maxextent = PAR_MAX(aabb[2] - aabb[0], aabb[3] - aabb[1]);
    PARFLT padding = 0.05 * maxextent;
    FILE* svgfile = fopen(filename, "wt");
    fprintf(svgfile,
        "<svg viewBox='%f %f %f %f' width='700px' height='700px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g stroke-width='0.5' stroke-opacity='0.5' stroke='black' "
        "fill-opacity='0.2' fill='#2A8BB6'>\n"
        "<rect fill-opacity='0.1' stroke='none' fill='#2A8BB6' x='%f' y='%f' "
        "width='100%%' height='100%%'/>\n",
        aabb[0] - padding, aabb[1] - padding,
        aabb[2] - aabb[0] + 2 * padding, aabb[3] - aabb[1] + 2 * padding,
        aabb[0] - padding, aabb[1] - padding);
    PARFLT const* xyr = bubbles->xyr;
    for (PARINT i = 0; i < bubbles->count; i++, xyr += 3) {
        fprintf(svgfile, "<circle stroke-width='%f' cx='%f' cy='%f' r='%f'/>\n",
            xyr[2] * 0.01, xyr[0], xyr[1], xyr[2]);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
}

#undef PARINT
#undef PARFLT
#endif // PAR_BUBBLES_IMPLEMENTATION
