// BUBBLES :: https://github.com/prideout/par
// Simple C library for packing circles into hierarchical (or flat) diagrams.
//
// Implements "Visualization of Large Hierarchical Data by Circle Packing" from
// Wang et al (2006).
//
// Also contains an implementation of Emo Welzl's "Smallest enclosing disks"
// algorithm (1991) for enclosing points with circles, and a related method
// from Mike Bostock for enclosing circles with the smallest possible circle.
//
// The API is divided into four sections:
//
//   - Enclosing.  Compute the smallest bounding circle for points or circles.
//   - Packing.    Pack circles together, or into other circles.
//   - Queries.    Given a touch point, pick a circle from a hierarchy, etc.
//   - Deep Zoom.  Uses relative coordinate systems to support arbitrary depth.
//
// In addition to the comment block above each function declaration, the API
// has informal documentation here:
//
//     https://prideout.net/bubbles
//
// Distributed under the MIT License, see bottom of file.

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

// This must be "float" or "double" or "long double". Note that you should not
// need high precision if you use the relative coordinate systems API.
#ifndef PAR_BUBBLES_FLT
#define PAR_BUBBLES_FLT double
#endif

// Enclosing / Touching --------------------------------------------------------

// Read an array of (x,y) coordinates, write a single 3-tuple (x,y,radius).
void par_bubbles_enclose_points(PAR_BUBBLES_FLT const* xy, PAR_BUBBLES_INT npts,
    PAR_BUBBLES_FLT* result);

// Read an array of 3-tuples (x,y,radius), write a 3-tuple (x,y,radius).
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

// Returns the smallest circle that intersects the three specified circles.
// This is the problem of Problem of Apollonius!
void par_bubbles_touch_three_disks(PAR_BUBBLES_FLT const* xyr1,
    PAR_BUBBLES_FLT const* xyr2, PAR_BUBBLES_FLT const* xyr3,
    PAR_BUBBLES_FLT* result);

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

// Returns a pointer to a list of children nodes.
void par_bubbles_get_children(par_bubbles_t const* bubbles, PAR_BUBBLES_INT idx,
    PAR_BUBBLES_INT** pchildren, PAR_BUBBLES_INT* nchildren);

// Returns the given node's parent, or 0 if it's the root.
PAR_BUBBLES_INT par_bubbles_get_parent(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT idx);

// Finds the height of the tree and returns one of its deepest leaves.
void par_bubbles_get_maxdepth(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT* maxdepth, PAR_BUBBLES_INT* leaf);

// Finds the height of the tree at a certain node.
PAR_BUBBLES_INT par_bubbles_get_depth(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT node);

// Returns a 4-tuple (min xy, max xy) for the given node.
void par_bubbles_compute_aabb_for_node(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT node, PAR_BUBBLES_FLT* aabb);

// Find the deepest node that is an ancestor of both A and B.  Classic!
PAR_BUBBLES_INT par_bubbles_lowest_common_ancestor(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT node_a, PAR_BUBBLES_INT node_b);

// Deep Zoom -------------------------------------------------------------------

// Similar to hpack, but maintains precision by storing disk positions within
// the local coordinate system of their parent. After calling this function,
// clients can use cull_local to flatten the coordinate systems.
par_bubbles_t* par_bubbles_hpack_local(PAR_BUBBLES_INT* nodes,
    PAR_BUBBLES_INT nnodes);

// Similar to par_bubbles_cull, but takes a root node rather than an AABB,
// and returns a result within the local coordinate system of the new root.
// In other words, the new root will have radius 1, centered at (0,0).  The
// minradius is also expressed in this coordinate system.
par_bubbles_t* par_bubbles_cull_local(par_bubbles_t const* src,
    PAR_BUBBLES_FLT const* aabb, PAR_BUBBLES_FLT minradius,
    PAR_BUBBLES_INT root, par_bubbles_t* dst);

// Finds the smallest node in the given bubble diagram that completely encloses
// the given axis-aligned bounding box (min xy, max xy).  The AABB coordinates
// are expressed in the local coordinate system of the given root node.
PAR_BUBBLES_INT par_bubbles_find_local(par_bubbles_t const* src,
    PAR_BUBBLES_FLT const* aabb, PAR_BUBBLES_INT root);

// Similar to pick, but expects (x,y) to be in the coordinate system of the
// given root node.
PAR_BUBBLES_INT par_bubbles_pick_local(par_bubbles_t const*, PAR_BUBBLES_FLT x,
    PAR_BUBBLES_FLT y, PAR_BUBBLES_INT root, PAR_BUBBLES_FLT minradius);

// Obtains the scale and translation (which should be applied in that order)
// that can move a point from the node0 coord system to the node1 coord system.
// The "xform" argument should point to three floats, which will be populated
// with: x translation, y translation, and scale.
bool par_bubbles_transform_local(par_bubbles_t const* bubbles,
    PAR_BUBBLES_FLT* xform, PAR_BUBBLES_INT node0, PAR_BUBBLES_INT node1);

// Dump out a SVG file for diagnostic purposes.
void par_bubbles_export_local(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT idx, char const* filename);

typedef enum {
    PAR_BUBBLES_FILTER_DEFAULT,
    PAR_BUBBLES_FILTER_DISCARD_LAST_CHILD,
    PAR_BUBBLES_FILTER_KEEP_ONLY_LAST_CHILD
} par_bubbles_filter;

// Special-case function that affects the behavior of subsequent calls to
// cull_local.  Allows clients to filter the children list of each non-leaf
// node, which is especially useful when using placeholder bubbles for labels.
void par_bubbles_set_filter(par_bubbles_t* src, par_bubbles_filter f);

typedef enum {
    PAR_BUBBLES_HORIZONTAL,
    PAR_BUBBLES_VERTICAL
} par_bubbles_orientation;

// Sets some global state that affect subsequent calls to hpack. The first two
// children can either be placed horizontally (default) or vertically.  The
// effect of this is subtle, since overall layout is obviously circular.
void par_bubbles_set_orientation(par_bubbles_orientation );

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

#ifdef __cplusplus
}
#endif

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

static par_bubbles_orientation par_bubbles__ostate = PAR_BUBBLES_HORIZONTAL;

typedef struct {
    PARINT prev;
    PARINT next;
} par_bubbles__node;

typedef struct {
    PARFLT* xyr;                 // results array
    PARINT count;                // client-provided count
    PARINT* ids;                 // populated by par_bubbles_cull
    PARFLT const* radiuses;      // client-provided radius list
    par_bubbles__node* chain;    // counterclockwise enveloping chain
    PARINT const* graph_parents; // client-provided parent indices
    PARINT* graph_children;      // flat list of children indices
    PARINT* graph_heads;         // list of "pointers" to first child
    PARINT* graph_tails;         // list of "pointers" to one-past-last child
    PARINT npacked;
    PARINT maxwidth;
    PARINT capacity;
    par_bubbles_filter filter;
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
    PARFLT x0, y0, x1, y1;
    if (par_bubbles__ostate == PAR_BUBBLES_HORIZONTAL) {
        x0 = -radii[0];
        y0 = 0;
        x1 = radii[1];
        y1 = 0;
    } else {
        x0 = 0;
        y0 = -radii[0];
        x1 = 0;
        y1 = radii[1];
    }
    *xyr++ = x0;
    *xyr++ = y0;
    *xyr++ = *radii++;
    if (bubbles->count == ++bubbles->npacked) {
        return;
    }
    *xyr++ = x1;
    *xyr++ = y1;
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

static void par__disk_from_two_points(PARFLT const* xy1, PARFLT const* xy2,
    PARFLT* result)
{
    PARFLT dx = xy1[0] - xy2[0];
    PARFLT dy = xy1[1] - xy2[1];
    result[0] = 0.5 * (xy1[0] + xy2[0]);
    result[1] = 0.5 * (xy1[1] + xy2[1]);
    result[2] = sqrt(dx * dx + dy * dy) / 2.0;
}

static PARINT par__disk_contains_point(PARFLT const* xyr, PARFLT const* xy)
{
    PARFLT dx = xyr[0] - xy[0];
    PARFLT dy = xyr[1] - xy[1];
    return dx * dx + dy * dy <= PAR_SQR(xyr[2]);
}

static void par__easydisk_from_points(PARFLT* disk, PARFLT const* edgepts,
    PARINT nedgepts)
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
    par__disk_from_two_points(edgepts, edgepts + 2, disk);
    if (nedgepts == 2 || par__disk_contains_point(disk, edgepts + 4)) {
        return;
    }
    par__disk_from_two_points(edgepts, edgepts + 4, disk);
    if (par__disk_contains_point(disk, edgepts + 2)) {
        return;
    }
    par__disk_from_two_points(edgepts + 2, edgepts + 4, disk);
    if (par__disk_contains_point(disk, edgepts)) {
        return;
    }
    par_bubbles_touch_three_points(edgepts, disk);
}

static void par__minidisk_points(PARFLT* disk, PARFLT const* pts, PARINT npts,
    PARFLT const* edgepts, PARINT nedgepts)
{
    if (npts == 0 || nedgepts == 3) {
        par__easydisk_from_points(disk, edgepts, nedgepts);
        return;
    }
    PARFLT const* pt = pts + (--npts) * 2;
    par__minidisk_points(disk, pts, npts, edgepts, nedgepts);
    if (!par__disk_contains_point(disk, pt)) {
        PARFLT edgepts1[6];
        for (PARINT i = 0; i < nedgepts * 2; i += 2) {
            edgepts1[i] = edgepts[i];
            edgepts1[i + 1] = edgepts[i + 1];
        }
        edgepts1[2 * nedgepts] = pt[0];
        edgepts1[2 * nedgepts + 1] = pt[1];
        par__minidisk_points(disk, pts, npts, edgepts1, ++nedgepts);
    }
}

// Returns true if the specified circle1 contains the specified circle2.
static bool par__disk_contains_disk(PARFLT const* xyr1, PARFLT const* xyr2)
{
    PARFLT xc0 = xyr1[0] - xyr2[0];
    PARFLT yc0 = xyr1[1] - xyr2[1];
    return sqrt(xc0 * xc0 + yc0 * yc0) < xyr1[2] - xyr2[2] + 1e-6;
}

// Returns the smallest circle that intersects the two specified circles.
static void par__disk_from_two_disks(PARFLT const* xyr1, PARFLT const* xyr2,
    PARFLT* result)
{
    PARFLT x1 = xyr1[0], y1 = xyr1[1], r1 = xyr1[2];
    PARFLT x2 = xyr2[0], y2 = xyr2[1], r2 = xyr2[2];
    PARFLT x12 = x2 - x1, y12 = y2 - y1, r12 = r2 - r1;
    PARFLT l = sqrt(x12 * x12 + y12 * y12);
    result[0] = (x1 + x2 + x12 / l * r12) / 2;
    result[1] = (y1 + y2 + y12 / l * r12) / 2;
    result[2] = (l + r1 + r2) / 2;
}

static void par__easydisk_from_disks(PARFLT* disk, PARFLT const* edgedisks,
    PARINT nedgedisks)
{
    assert(nedgedisks <= 3);
    if (nedgedisks == 0) {
        disk[0] = 0;
        disk[1] = 0;
        disk[2] = 0;
        return;
    }
    if (nedgedisks == 1) {
        disk[0] = edgedisks[0];
        disk[1] = edgedisks[1];
        disk[2] = edgedisks[2];
        return;
    }
    if (nedgedisks == 2) {
        par__disk_from_two_disks(edgedisks, edgedisks + 3, disk);
        return;
    }
    par_bubbles_touch_three_disks(edgedisks, edgedisks + 3, edgedisks + 6,
        disk);
}

static void par__minidisk_disks(PARFLT* result, PARFLT const* disks,
    PARINT ndisks, PARFLT const* edgedisks, PARINT nedgedisks)
{
    if (ndisks == 0 || nedgedisks == 3) {
        par__easydisk_from_disks(result, edgedisks, nedgedisks);
        return;
    }
    PARFLT const* disk = disks + (--ndisks) * 3;
    par__minidisk_disks(result, disks, ndisks, edgedisks, nedgedisks);
    if (!par__disk_contains_disk(result, disk)) {
        PARFLT edgedisks1[9];
        for (PARINT i = 0; i < nedgedisks * 3; i += 3) {
            edgedisks1[i] = edgedisks[i];
            edgedisks1[i + 1] = edgedisks[i + 1];
            edgedisks1[i + 2] = edgedisks[i + 2];
        }
        edgedisks1[3 * nedgedisks] = disk[0];
        edgedisks1[3 * nedgedisks + 1] = disk[1];
        edgedisks1[3 * nedgedisks + 2] = disk[2];
        par__minidisk_disks(result, disks, ndisks, edgedisks1, ++nedgedisks);
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
    par__minidisk_points(result, xy, npts, 0, 0);
}

void par_bubbles_enclose_disks(PARFLT const* xyr, PARINT ndisks, PARFLT* result)
{
    par__minidisk_disks(result, xyr, ndisks, 0, 0);
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

void par_bubbles_touch_three_disks(PARFLT const* xyr1, PARFLT const* xyr2,
    PARFLT const* xyr3, PARFLT* result)
{
    PARFLT x1 = xyr1[0], y1 = xyr1[1], r1 = xyr1[2],
        x2 = xyr2[0], y2 = xyr2[1], r2 = xyr2[2],
        x3 = xyr3[0], y3 = xyr3[1], r3 = xyr3[2],
        a2 = 2 * (x1 - x2),
        b2 = 2 * (y1 - y2),
        c2 = 2 * (r2 - r1),
        d2 = x1 * x1 + y1 * y1 - r1 * r1 - x2 * x2 - y2 * y2 + r2 * r2,
        a3 = 2 * (x1 - x3),
        b3 = 2 * (y1 - y3),
        c3 = 2 * (r3 - r1),
        d3 = x1 * x1 + y1 * y1 - r1 * r1 - x3 * x3 - y3 * y3 + r3 * r3,
        ab = a3 * b2 - a2 * b3,
        xa = (b2 * d3 - b3 * d2) / ab - x1,
        xb = (b3 * c2 - b2 * c3) / ab,
        ya = (a3 * d2 - a2 * d3) / ab - y1,
        yb = (a2 * c3 - a3 * c2) / ab,
        A = xb * xb + yb * yb - 1,
        B = 2 * (xa * xb + ya * yb + r1),
        C = xa * xa + ya * ya - r1 * r1,
        r = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
    result[0] = xa + xb * r + x1;
    result[1] = ya + yb * r + y1;
    result[2] = r;
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

// Assigns a radius to every node according to its number of descendants.
void par_bubbles__generate_radii(par_bubbles__t* bubbles,
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
        par_bubbles__generate_radii(bubbles, worker, child);
        bubbles->xyr[pr] += bubbles->xyr[child * 3 + 2];
    }
    // The following square root seems to produce a nicer, more space-filling,
    // distribution of radiuses  in randomly-generated trees.
    bubbles->xyr[pr] = sqrtf(bubbles->xyr[pr]);
}

void par_bubbles__hpack(par_bubbles__t* bubbles, par_bubbles__t* worker,
    PARINT parent, bool local)
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
    const PARFLT PAR_HPACK_PADDING1 = 0.15;
    const PARFLT PAR_HPACK_PADDING2 = 0.025;
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

    // Transform the children to fit nicely into either (a) the unit circle,
    // or (b) their parent.  The former is used if "local" is true.
    PARFLT scale, tx, ty;
    if (local) {
        scale = 1.0 / cr;
        tx = 0;
        ty = 0;
    } else {
        scale = pr / cr;
        tx = px;
        ty = py;
    }
    PARFLT const* src = worker->xyr;
    for (PARINT cindex = head; cindex != tail; cindex++, src += 3) {
        PARFLT* dst = bubbles->xyr + 3 * bubbles->graph_children[cindex];
        dst[0] = tx + scale * (src[0] - cx);
        dst[1] = ty + scale * (src[1] - cy);
        dst[2] = scale * (src[2] - scaled_padding);
    }

    // Recursion.  TODO: It might be better to use our own stack here.
    for (PARINT cindex = head; cindex != tail; cindex++) {
        par_bubbles__hpack(bubbles, worker, bubbles->graph_children[cindex],
            local);
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
        par_bubbles__generate_radii(bubbles, worker, 0);
        bubbles->xyr[0] = 0;
        bubbles->xyr[1] = 0;
        bubbles->xyr[2] = radius;
        par_bubbles__hpack(bubbles, worker, 0, false);
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
        "<svg viewBox='%f %f %f %f' width='640px' height='640px' "
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

PARINT par_bubbles_get_parent(par_bubbles_t const* pbubbles, PARINT node)
{
    par_bubbles__t const* bubbles = (par_bubbles__t const*) pbubbles;
    return bubbles->graph_parents[node];
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

PARINT par_bubbles_lowest_common_ancestor(par_bubbles_t const* bubbles,
    PARINT node_a, PARINT node_b)
{
    if (node_a == node_b) {
        return node_a;
    }
    par_bubbles__t const* src = (par_bubbles__t const*) bubbles;
    PARINT depth_a = par_bubbles_get_depth(bubbles, node_a);
    PARINT* chain_a = PAR_MALLOC(PARINT, depth_a);
    for (PARINT i = depth_a - 1; i >= 0; i--) {
        chain_a[i] = node_a;
        node_a = src->graph_parents[node_a];
    }
    PARINT depth_b = par_bubbles_get_depth(bubbles, node_b);
    PARINT* chain_b = PAR_MALLOC(PARINT, depth_b);
    for (PARINT i = depth_b - 1; i >= 0; i--) {
        chain_b[i] = node_b;
        node_b = src->graph_parents[node_b];
    }
    PARINT lca = 0;
    for (PARINT i = 1; i < PAR_MIN(depth_a, depth_b); i++) {
        if (chain_a[i] != chain_b[i]) {
            break;
        }
        lca = chain_a[i];
    }
    PAR_FREE(chain_a);
    PAR_FREE(chain_b);
    return lca;
}

void par_bubbles_export_local(par_bubbles_t const* bubbles,
    PAR_BUBBLES_INT root, char const* filename)
{
    par_bubbles_t* clone = par_bubbles_cull_local(bubbles, 0, 0, root, 0);
    FILE* svgfile = fopen(filename, "wt");
    fprintf(svgfile,
        "<svg viewBox='%f %f %f %f' width='640px' height='640px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g stroke-width='0.5' stroke-opacity='0.5' stroke='black' "
        "fill-opacity='0.2' fill='#2A8BB6'>\n"
        "<rect fill-opacity='0.1' stroke='none' fill='#2AB68B' x='%f' y='%f' "
        "width='100%%' height='100%%'/>\n",
        -1.0, -1.0, 2.0, 2.0, -1.0, -1.0);
    PARFLT const* xyr = clone->xyr;
    for (PARINT i = 0; i < clone->count; i++, xyr += 3) {
        fprintf(svgfile, "<circle stroke-width='%f' cx='%f' cy='%f' r='%f'/>\n",
            xyr[2] * 0.01, xyr[0], xyr[1], xyr[2]);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
    par_bubbles_free_result(clone);
}

void par_bubbles_set_filter(par_bubbles_t* bubbles, par_bubbles_filter f)
{
    par_bubbles__t* src = (par_bubbles__t*) bubbles;
    src->filter = f;
}

static void par_bubbles__copy_disk_local(par_bubbles__t const* src,
    par_bubbles__t* dst, PARINT parent, PARFLT const* xform)
{
    PARINT i = dst->count++;
    if (dst->capacity < dst->count) {
        dst->capacity = PAR_MAX(16, dst->capacity) * 2;
        dst->xyr = PAR_REALLOC(PARFLT, dst->xyr, 3 * dst->capacity);
        dst->ids = PAR_REALLOC(PARINT, dst->ids, dst->capacity);
    }
    PARFLT const* xyr = src->xyr + parent * 3;
    dst->xyr[i * 3] = xyr[0] * xform[2] + xform[0];
    dst->xyr[i * 3 + 1] = xyr[1] * xform[2] + xform[1];
    dst->xyr[i * 3 + 2] = xyr[2] * xform[2];
    dst->ids[i] = parent;
}

static void par_bubbles__cull_local(par_bubbles__t const* src,
    PARFLT const* aabb, PARFLT const* xform, PARFLT minradius,
    par_bubbles__t* dst, PARINT parent)
{
    PARFLT const* xyr = src->xyr + parent * 3;
    PARFLT child_xform[3] = {
        xform[0] + xform[2] * xyr[0],
        xform[1] + xform[2] * xyr[1],
        xform[2] * xyr[2]
    };
    if (aabb && !par_bubbles_check_aabb(child_xform, aabb)) {
        return;
    }
    if (child_xform[2] < minradius) {
        return;
    }
    par_bubbles__copy_disk_local(src, dst, parent, xform);
    xform = child_xform;
    PARINT head = src->graph_heads[parent];
    PARINT tail = src->graph_tails[parent];
    if (src->filter == PAR_BUBBLES_FILTER_DISCARD_LAST_CHILD) {
        tail--;
    } else if (src->filter == PAR_BUBBLES_FILTER_KEEP_ONLY_LAST_CHILD) {
        head = PAR_MAX(head, tail - 1);
    }
    for (PARINT cindex = head; cindex < tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        par_bubbles__cull_local(src, aabb, xform, minradius, dst, child);
    }
}

par_bubbles_t* par_bubbles_cull_local(par_bubbles_t const* psrc,
    PAR_BUBBLES_FLT const* aabb, PAR_BUBBLES_FLT minradius,
    PAR_BUBBLES_INT root, par_bubbles_t* pdst)
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
    PARFLT xform[3] = {0, 0, 1};
    par_bubbles__copy_disk_local(src, dst, root, xform);
    dst->xyr[0] = dst->xyr[1] = 0;
    dst->xyr[2] = 1;
    PARINT head = src->graph_heads[root];
    PARINT tail = src->graph_tails[root];
    if (src->filter == PAR_BUBBLES_FILTER_DISCARD_LAST_CHILD) {
        tail--;
    } else if (src->filter == PAR_BUBBLES_FILTER_KEEP_ONLY_LAST_CHILD) {
        head = PAR_MAX(head, tail - 1);
    }
    for (PARINT cindex = head; cindex < tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        par_bubbles__cull_local(src, aabb, xform, minradius, dst, child);
    }
    return pdst;
}

par_bubbles_t* par_bubbles_hpack_local(PARINT* nodes, PARINT nnodes)
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
        par_bubbles__generate_radii(bubbles, worker, 0);
        bubbles->xyr[0] = 0;
        bubbles->xyr[1] = 0;
        bubbles->xyr[2] = 1;
        par_bubbles__hpack(bubbles, worker, 0, true);
        par_bubbles_free_result((par_bubbles_t*) worker);
    }
    return (par_bubbles_t*) bubbles;
}

static bool par_bubbles__disk_encloses_aabb(PAR_BUBBLES_FLT cx,
    PAR_BUBBLES_FLT cy, PAR_BUBBLES_FLT r, PAR_BUBBLES_FLT const* aabb)
{
    PAR_BUBBLES_FLT x, y;
    PAR_BUBBLES_FLT r2 = r * r;
    x = aabb[0]; y = aabb[1];
    if (PAR_SQR(x - cx) + PAR_SQR(y - cy) > r2) {
        return false;
    }
    x = aabb[2]; y = aabb[1];
    if (PAR_SQR(x - cx) + PAR_SQR(y - cy) > r2) {
        return false;
    }
    x = aabb[0]; y = aabb[3];
    if (PAR_SQR(x - cx) + PAR_SQR(y - cy) > r2) {
        return false;
    }
    x = aabb[2]; y = aabb[3];
    return PAR_SQR(x - cx) + PAR_SQR(y - cy) <= r2;
}

static bool par_bubbles__get_local(par_bubbles__t const* src, PARFLT* xform,
    PARINT parent, PARINT node);

static bool par_bubbles_transform_parent(par_bubbles__t const* src,
    PARFLT* xform, PARINT node0)
{
    PARINT node1 = src->graph_parents[node0];
    xform[0] = 0;
    xform[1] = 0;
    xform[2] = 1;

    PARINT head = src->graph_heads[node1];
    PARINT tail = src->graph_tails[node1];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        if (par_bubbles__get_local(src, xform, child, node0)) {
            return true;
        }
    }
    return false;
}

PARINT par_bubbles__find_local(par_bubbles__t const* src,
    PARFLT const* xform, PARFLT const* aabb, PARINT parent)
{
    PARFLT const* xyr = src->xyr + parent * 3;
    PARFLT child_xform[3] = {
        xform[2] * xyr[0] + xform[0],
        xform[2] * xyr[1] + xform[1],
        xform[2] * xyr[2]
    };
    xform = child_xform;
    if (!par_bubbles__disk_encloses_aabb(xform[0], xform[1], xform[2], aabb)) {
        return -1;
    }
    PARFLT maxrad = 0;
    PARINT head = src->graph_heads[parent];
    PARINT tail = src->graph_tails[parent];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        PARFLT const* xyr = src->xyr + child * 3;
        maxrad = PAR_MAX(maxrad, xyr[2]);
    }
    PARFLT maxext = PAR_MAX(aabb[2] - aabb[0], aabb[3] - aabb[1]);
    if (2 * maxrad < maxext) {
        return parent;
    }
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        PARINT cresult = par_bubbles__find_local(src, xform, aabb, child);
        if (cresult > -1) {
            return cresult;
        }
    }
    return parent;
}

// This finds the deepest node that completely encloses the box.
PARINT par_bubbles_find_local(par_bubbles_t const* bubbles, PARFLT const* aabb,
    PARINT root)
{
    par_bubbles__t const* src = (par_bubbles__t const*) bubbles;

    // Since the aabb is expressed in the coordinate system of the given root,
    // we can do a trivial rejection right away, using the unit circle.
    if (!par_bubbles__disk_encloses_aabb(0, 0, 1, aabb)) {
        if (root == 0) {
            return -1;
        }
        PARFLT xform[3];
        par_bubbles_transform_parent(src, xform, root);
        PARFLT width = aabb[2] - aabb[0];
        PARFLT height = aabb[3] - aabb[1];
        PARFLT cx = 0.5 * (aabb[0] + aabb[2]);
        PARFLT cy = 0.5 * (aabb[1] + aabb[3]);
        width *= xform[2];
        height *= xform[2];
        cx = cx * xform[2] + xform[0];
        cy = cy * xform[2] + xform[1];
        PARFLT new_aabb[4] = {
            cx - width * 0.5,
            cy - height * 0.5,
            cx + width * 0.5,
            cy + height * 0.5
        };
        PARINT parent = src->graph_parents[root];
        return par_bubbles_find_local(bubbles, new_aabb, parent);
    }

    PARFLT maxrad = 0;
    PARINT head = src->graph_heads[root];
    PARINT tail = src->graph_tails[root];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        PARFLT const* xyr = src->xyr + child * 3;
        maxrad = PAR_MAX(maxrad, xyr[2]);
    }
    PARFLT maxext = PAR_MAX(aabb[2] - aabb[0], aabb[3] - aabb[1]);
    if (2 * maxrad < maxext) {
        return root;
    }

    PARFLT xform[3] = {0, 0, 1};
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        PARINT cresult = par_bubbles__find_local(src, xform, aabb, child);
        if (cresult > -1) {
            return cresult;
        }
    }
    return root;
}

// This could be implemented much more efficiently, but for now it simply
// calls find_local with a zero-size AABB, then ensures that the result
// has a radius that is greater than or equal to minradius.
PARINT par_bubbles_pick_local(par_bubbles_t const* bubbles, PARFLT x, PARFLT y,
    PARINT root, PARFLT minradius)
{
    par_bubbles__t const* src = (par_bubbles__t const*) bubbles;
    PARFLT aabb[] = { x, y, x, y };
    PARINT result = par_bubbles_find_local(bubbles, aabb, root);
    if (result == -1) {
        return result;
    }
    PARINT depth = par_bubbles_get_depth(bubbles, result);
    PARINT* chain = PAR_MALLOC(PARINT, depth);
    PARINT node = result;
    for (PARINT i = depth - 1; i >= 0; i--) {
        chain[i] = node;
        node = src->graph_parents[node];
    }
    PARFLT radius = 1;
    for (PARINT i = 1; i < depth; i++) {
        PARINT node = chain[i];
        radius *= src->xyr[node * 3 + 2];
        if (radius < minradius) {
            result = chain[i - 1];
            break;
        }
    }
    PAR_FREE(chain);
    return result;
}

static bool par_bubbles__get_local(par_bubbles__t const* src, PARFLT* xform,
    PARINT parent, PARINT node)
{
    PARFLT const* xyr = src->xyr + parent * 3;
    PARFLT child_xform[3] = {
        xform[2] * xyr[0] + xform[0],
        xform[2] * xyr[1] + xform[1],
        xform[2] * xyr[2]
    };
    if (parent == node) {
        xform[0] = child_xform[0];
        xform[1] = child_xform[1];
        xform[2] = child_xform[2];
        return true;
    }
    PARINT head = src->graph_heads[parent];
    PARINT tail = src->graph_tails[parent];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        if (par_bubbles__get_local(src, child_xform, child, node)) {
            xform[0] = child_xform[0];
            xform[1] = child_xform[1];
            xform[2] = child_xform[2];
            return true;
        }
    }
    return false;
}

// Obtains the scale and translation (which should be applied in that order)
// that can move a point from the node0 coord system to the node1 coord system.
// The "xform" argument should point to three floats, which will be populated
// with: x translation, y translation, and scale.
bool par_bubbles_transform_local(par_bubbles_t const* bubbles, PARFLT* xform,
    PARINT node0, PARINT node1)
{
    par_bubbles__t const* src = (par_bubbles__t const*) bubbles;
    xform[0] = 0;
    xform[1] = 0;
    xform[2] = 1;
    if (node0 == node1) {
        return true;
    }
    if (node1 == src->graph_parents[node0]) {
        return par_bubbles_transform_parent(src, xform, node0);
    }

    // First try the case where node1 is a descendant of node0
    PARINT head = src->graph_heads[node0];
    PARINT tail = src->graph_tails[node0];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        if (par_bubbles__get_local(src, xform, child, node1)) {
            float tx = xform[0];
            float ty = xform[1];
            float s = xform[2];
            xform[0] = -tx / s;
            xform[1] = -ty / s;
            xform[2] = 1.0 / s;
            return true;
        }
    }

    // Next, try the case where node0 is a descendant of node1
    head = src->graph_heads[node1];
    tail = src->graph_tails[node1];
    for (PARINT cindex = head; cindex != tail; cindex++) {
        PARINT child = src->graph_children[cindex];
        if (par_bubbles__get_local(src, xform, child, node0)) {
            return true;
        }
    }

    // If we reach here, then node0 is neither an ancestor nor a descendant, so
    // do something hacky and return false.  It would be best to find the lowest
    // common ancestor, but let's just assume the lowest common ancestor is 0.
    PARFLT xform2[3] = {0, 0, 1};
    par_bubbles_transform_local(bubbles, xform, node0, 0);
    par_bubbles_transform_local(bubbles, xform2, 0, node1);
    xform[0] *= xform2[2];
    xform[1] *= xform2[2];
    xform[2] *= xform2[2];
    xform[0] += xform2[0];
    xform[1] += xform2[1];

    return false;
}

void par_bubbles_set_orientation(par_bubbles_orientation ostate)
{
    par_bubbles__ostate = ostate;
}

#undef PARINT
#undef PARFLT
#endif // PAR_BUBBLES_IMPLEMENTATION
#endif // PAR_BUBBLES_H

// par_bubbles is distributed under the MIT license:
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
