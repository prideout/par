// SPRUNE :: https://github.com/prideout/par
// Sweep and Prune library for detecting axis-aligned box collisions in 2D.
//
// For an emscripten demo of this library, take a look at the following link.
//
//     https://prideout.net/d3cpp
//
// The axis-aligned bounding boxes are specified by (minx, miny, maxx, maxy).
// Simple usage example:
//
//     float boxes[] = {
//         0.10, 0.10, 0.30, 0.30, // box 0
//         0.20, 0.20, 0.40, 0.40, // box 1
//         0.60, 0.15, 0.70, 0.25, // box 2
//     };
//     int nboxes = 3;
//     par_sprune_context* ctx = par_sprune_overlap(boxes, nboxes, 0);
//     int const* pairs = ctx->collision_pairs;
//     for (int i = 0; i < ctx->ncollision_pairs * 2; i += 2) {
//         printf("box %d overlaps box %d\n", pairs[i], pairs[i + 1]);
//     }
//     par_sprune_free_context(ctx);
//
// Distributed under the MIT License, see bottom of file.

#ifndef PAR_SPRUNE_H
#define PAR_SPRUNE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>

#ifndef PAR_SPRUNE_INT
#define PAR_SPRUNE_INT int32_t
#endif

#ifndef PAR_SPRUNE_FLT
#define PAR_SPRUNE_FLT float
#endif

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

typedef struct {
    PAR_SPRUNE_INT const* const collision_pairs; // list of two-tuples
    PAR_SPRUNE_INT const ncollision_pairs;       // number of two-tuples
    PAR_SPRUNE_INT const* const culled;          // filled by par_sprune_cull
    PAR_SPRUNE_INT const nculled;                // set by par_sprune_cull
} par_sprune_context;

void par_sprune_free_context(par_sprune_context* context);

// Takes an array of 4-tuples (minx miny maxx maxy) and performs SaP. Populates
// "collision_pairs" and "ncollision_pairs".  Optionally takes an existing
// context to avoid memory churn; pass NULL for initial construction.
par_sprune_context* par_sprune_overlap(PAR_SPRUNE_FLT const* aabbs,
    PAR_SPRUNE_INT naabbs, par_sprune_context* previous);

// Reads new aabb data from the same pointer that was passed to the overlap
// function and refreshes the two relevant fields.  This function should
// only be used when the number of aabbs remains constant. If this returns
// false, no changes to the collision set were detected.
bool par_sprune_update(par_sprune_context* ctx);

// Examines all collision groups and creates a culling set such that no boxes
// would overlap if the culled boxes are removed.  When two boxes collide, the
// box that occurs earlier in the list is more likely to be culled. Populates
// the "culled" and "nculled" fields in par_sprune_context. This is useful for
// hiding labels in GIS applications.
void par_sprune_cull(par_sprune_context* context);

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#ifdef PAR_SPRUNE_IMPLEMENTATION
#define PARINT PAR_SPRUNE_INT
#define PARFLT PAR_SPRUNE_FLT

#include <stdlib.h>
#include <assert.h>

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

#ifndef PAR_ARRAY
#define PAR_ARRAY
#define pa_free(a) ((a) ? PAR_FREE(pa___raw(a)), 0 : 0)
#define pa_push(a, v) (pa___maybegrow(a, (int) 1), (a)[pa___n(a)++] = (v))
#define pa_count(a) ((a) ? pa___n(a) : 0)
#define pa_add(a, n) (pa___maybegrow(a, (int) n), pa___n(a) += (n))
#define pa_last(a) ((a)[pa___n(a) - 1])
#define pa_end(a) (a + pa_count(a))
#define pa_clear(arr) if (arr) pa___n(arr) = 0
#define pa___raw(a) ((int*) (a) -2)
#define pa___m(a) pa___raw(a)[0]
#define pa___n(a) pa___raw(a)[1]
#define pa___needgrow(a, n) ((a) == 0 || pa___n(a) + ((int) n) >= pa___m(a))
#define pa___maybegrow(a, n) (pa___needgrow(a, (n)) ? pa___grow(a, n) : 0)
#define pa___grow(a, n) (*((void**)& (a)) = pa___growf((void*) (a), (n), \
        sizeof(*(a))))

// ptr[-2] is capacity, ptr[-1] is size.
static void* pa___growf(void* arr, int increment, int itemsize)
{
    int dbl_cur = arr ? 2 * pa___m(arr) : 0;
    int min_needed = pa_count(arr) + increment;
    int m = dbl_cur > min_needed ? dbl_cur : min_needed;
    int* p = (int *) PAR_REALLOC(uint8_t, arr ? pa___raw(arr) : 0,
        itemsize * m + sizeof(int) * 2);
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

typedef struct {

    // Public:
    PARINT* collision_pairs;
    PARINT ncollision_pairs;
    PARINT* culled;
    PARINT nculled;

    // Private:
    PARFLT const* aabbs;
    PARINT naabbs;
    PARINT* sorted_indices[2];
    PARINT* pairs[2];

} par_sprune__context;

static inline int par_qsort_cmpswap(char *__restrict a, char *__restrict b,
    size_t w,
    int (*compar)(const void *_a, const void *_b,
    void *_arg),
    void *arg)
{
    char tmp, *end = a+w;
    if (compar(a, b, arg) > 0) {
        for(; a < end; a++, b++) { tmp = *a; *a = *b; *b = tmp; }
        return 1;
    }
    return 0;
}

// qsort doesn't take a context, so we have our own portable implementation.
// Parameters:
//     base is the array to be sorted
//     nel is the number of elements in the array
//     w is the size in bytes of each element of the array
//     compar is the comparison function
//     arg is a pointer to be passed to the comparison function
//
static inline void par_qsort(
    void *base,
    size_t nel,
    size_t w,
    int (*compar)(const void *_a, const void *_b, void *_arg),
    void *arg)
{
    char *b = (char*) base, *end = (char*) (b + nel * w);
    if (nel < 7) {
        char *pi, *pj;
        for (pi = b+w; pi < end; pi += w) {
            for (pj = pi; pj > b && par_qsort_cmpswap(pj-w, pj, w, compar, arg);
                pj -= w) {}
        }
        return;
    }
    char *x, *y, *xend, ch;
    char *pl, *pr;
    char *last = b+w*(nel-1), *tmp;
    char *l[3];
    l[0] = b;
    l[1] = b+w*(nel/2);
    l[2] = last;
    if (compar(l[0],l[1],arg) > 0) {
        tmp=l[0]; l[0]=l[1]; l[1]=tmp;
    }
    if (compar(l[1],l[2],arg) > 0) {
        tmp=l[1]; l[1]=l[2]; l[2]=tmp;
        if (compar(l[0],l[1],arg) > 0) {
            tmp=l[0]; l[0]=l[1]; l[1]=tmp;
        }
    }
    for(x = l[1], y = last, xend = x+w; x<xend; x++, y++) {
        ch = *x; *x = *y; *y = ch;
    }
    pl = b;
    pr = last;
    while (pl < pr) {
        for (; pl < pr; pl += w) {
            if (par_qsort_cmpswap(pl, pr, w, compar, arg)) {
                pr -= w;
                break;
            }
        }
        for (; pl < pr; pr -= w) {
            if (par_qsort_cmpswap(pl, pr, w, compar, arg)) {
                pl += w;
                break;
            }
        }
    }
    par_qsort(b, (pl-b) / w, w, compar, arg);
    par_qsort(pl+w, (end - (pl+w)) / w, w, compar, arg);
}

void par_sprune_free_context(par_sprune_context* context)
{
    par_sprune__context* ctx = (par_sprune__context*) context;
    pa_free(ctx->sorted_indices[0]);
    pa_free(ctx->sorted_indices[1]);
    pa_free(ctx->pairs[0]);
    pa_free(ctx->pairs[1]);
    pa_free(ctx->collision_pairs);
    PAR_FREE(ctx);
}

static void par_sprune__remove(PARINT* arr, PARINT val)
{
    int i = pa_count(arr) - 1;
    for (; i >= 0; i--) {
        if (arr[i] == val) {
            break;
        }
    }
    assert(i >= 0);
    for (++i; i < pa_count(arr); i++) {
        PAR_SWAP(PARINT, arr[i - 1], arr[i]);
    }
    pa___n(arr)--;
}

typedef struct {
    PARFLT const* aabbs;
} par__sprune_sorter;

static int par__cmpinds(const void* pa, const void* pb, void* psorter)
{
    PARINT a = *((const PARINT*) pa);
    PARINT b = *((const PARINT*) pb);
    par__sprune_sorter* sorter = (par__sprune_sorter*) psorter;
    PARFLT const* aabbs = sorter->aabbs;
    PARFLT vala = aabbs[a];
    PARFLT valb = aabbs[b];
    if (vala > valb) return 1;
    if (vala < valb) return -1;
    if (a > b) return 1;
    if (a < b) return -1;
    return 0;
}

static int par__cmppairs(const void* pa, const void* pb, void* unused)
{
    PARINT a = *((const PARINT*) pa);
    PARINT b = *((const PARINT*) pb);
    if (a > b) return 1;
    if (a < b) return -1;
    a = *(1 + (const PARINT*) pa);
    b = *(1 + (const PARINT*) pb);
    if (a > b) return 1;
    if (a < b) return -1;
    return 0;
}

static int par__cmpfind(const void* pa, const void* pb)
{
    PARINT a = *((const PARINT*) pa);
    PARINT b = *((const PARINT*) pb);
    if (a > b) return 1;
    if (a < b) return -1;
    a = *(1 + (const PARINT*) pa);
    b = *(1 + (const PARINT*) pb);
    if (a > b) return 1;
    if (a < b) return -1;
    return 0;
}

par_sprune_context* par_sprune_overlap(PARFLT const* aabbs, PARINT naabbs,
    par_sprune_context* previous)
{
    par_sprune__context* ctx = (par_sprune__context*) previous;
    if (!ctx) {
        ctx = PAR_CALLOC(par_sprune__context, 1);
    }
    ctx->aabbs = aabbs;
    ctx->naabbs = naabbs;
    for (int axis = 0; axis < 2; axis++) {
        pa_clear(ctx->sorted_indices[axis]);
        pa_add(ctx->sorted_indices[axis], naabbs * 2);
        pa_clear(ctx->pairs[axis]);
    }
    for (PARINT i = 0; i < naabbs; i++) {
        ctx->sorted_indices[0][i * 2 + 0] = i * 4 + 0;
        ctx->sorted_indices[1][i * 2 + 0] = i * 4 + 1;
        ctx->sorted_indices[0][i * 2 + 1] = i * 4 + 2;
        ctx->sorted_indices[1][i * 2 + 1] = i * 4 + 3;
    }
    par__sprune_sorter sorter;
    sorter.aabbs = ctx->aabbs;
    PARINT* active = 0;

    // Sweep a plane first across the X-axis, then down through the Y-axis.

    for (int axis = 0; axis < 2; axis++) {
        PARINT** pairs = &ctx->pairs[axis];
        PARINT* indices = ctx->sorted_indices[axis];
        par_qsort(indices, naabbs * 2, sizeof(PARINT), par__cmpinds, &sorter);
        pa_clear(active);
        for (PARINT i = 0; i < naabbs * 2; i++) {
            PARINT fltindex = indices[i];
            PARINT boxindex = fltindex / 4;
            bool ismin = ((fltindex - axis) % 4) == 0;
            if (ismin) {
                for (int j = 0; j < pa_count(active); j++) {
                    pa_push(*pairs, active[j]);
                    pa_push(*pairs, boxindex);
                    pa_push(*pairs, boxindex);
                    pa_push(*pairs, active[j]);
                }
                pa_push(active, boxindex);
            } else {
                par_sprune__remove(active, boxindex);
            }
        }
    }

    // Sort the Y-axis collision pairs to make it easier to intersect it
    // with the set of X-axis collision pairs.  We also sort the X-axis
    // pairs because it's required for subsequent calls to par_sprune_update.

    PARINT* xpairs = ctx->pairs[0];
    PARINT* ypairs = ctx->pairs[1];
    int nxpairs = pa_count(xpairs) / 2;
    int nypairs = pa_count(ypairs) / 2;
    int pairsize = 2 * sizeof(PARINT);
    pa_free(active);
    par_qsort(xpairs, nxpairs, pairsize, par__cmppairs, 0);
    par_qsort(ypairs, nypairs, pairsize, par__cmppairs, 0);
    pa_clear(ctx->collision_pairs);

    // Find the intersection of X-axis overlaps and Y-axis overlaps.

    for (int i = 0; i < pa_count(xpairs); i += 2) {
        PARINT* key = xpairs + i;
        if (key[1] < key[0]) {
            continue;
        }
        void* found = bsearch(key, ypairs, nypairs, pairsize, par__cmpfind);
        if (found) {
            pa_push(ctx->collision_pairs, key[0]);
            pa_push(ctx->collision_pairs, key[1]);
        }
    }
    ctx->ncollision_pairs = pa_count(ctx->collision_pairs) / 2;
    return (par_sprune_context*) ctx;
}

bool par_sprune_update(par_sprune_context* context)
{
    par_sprune__context* ctx = (par_sprune__context*) context;
    PARINT* collision_pairs = ctx->collision_pairs;
    PARINT ncollision_pairs = ctx->ncollision_pairs;
    ctx->collision_pairs = 0;
    par_sprune_overlap(ctx->aabbs, ctx->naabbs, context);
    bool dirty = ncollision_pairs != ctx->ncollision_pairs;
    if (!dirty) {
        int pairsize = 2 * sizeof(PARINT);
        for (int i = 0; i < ctx->ncollision_pairs; i += 2) {
            PARINT* key = ctx->collision_pairs + i;
            if (!bsearch(key, collision_pairs, ncollision_pairs,
                pairsize, par__cmpfind)) {
                dirty = true;
                break;
            }
        }
    }
    pa_free(collision_pairs);
    return dirty;
}

bool par_sprune__is_culled(par_sprune__context* ctx, PARINT key)
{
    for (int i = 0; i < pa_count(ctx->culled); i++) {
        if (key == ctx->culled[i]) {
            return true;
        }
    }
    return false;
}

static int par__cmpfindsingle(const void* pa, const void* pb)
{
    PARINT a = *((const PARINT*) pa);
    PARINT b = *((const PARINT*) pb);
    if (a > b) return 1;
    if (a < b) return -1;
    return 0;
}

void par_sprune_cull(par_sprune_context* context)
{
    par_sprune__context* ctx = (par_sprune__context*) context;
    pa_clear(ctx->culled);
    PARINT* collision_pairs = ctx->collision_pairs;
    PARINT ncollision_pairs = ctx->ncollision_pairs;
    int pairsize = 2 * sizeof(PARINT);
    for (int i = 0; i < ctx->naabbs; i++) {
        PARINT* found = (PARINT*) bsearch(&i, collision_pairs, ncollision_pairs,
            pairsize, par__cmpfindsingle);
        if (!found) {
            continue;
        }
        if (!par_sprune__is_culled(ctx, found[0]) &&
            !par_sprune__is_culled(ctx, found[1])) {
            pa_push(ctx->culled, found[0]);
        }
    }
    ctx->nculled = pa_count(ctx->culled);
}

#undef PARINT
#undef PARFLT
#endif // PAR_SPRUNE_IMPLEMENTATION
#endif // PAR_SPRUNE_H

// par_sprune is distributed under the MIT license:
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
