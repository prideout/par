// BLUENOISE :: https://github.com/prideout/par
// Generator for infinite 2D point sequences using Recursive Wang Tiles.
//
// In addition to this source code, you'll need to download one of the following
// tilesets, the first being 2 MB while the other is 257 KB. The latter cheats
// by referencing the point sequence from the first tile for all 8 tiles. This
// obviously produces poor results, but in many contexts, it isn't noticeable.
//
//     https://prideout.net/assets/bluenoise.bin
//     https://prideout.net/assets/bluenoise.trimmed.bin
//
// The code herein is an implementation of the algorithm described in:
//
//     Recursive Wang Tiles for Real-Time Blue Noise
//     Johannes Kopf, Daniel Cohen-Or, Oliver Deussen, Dani Lischinski
//     ACM Transactions on Graphics 25, 3 (Proc. SIGGRAPH 2006)
//
// If you use this software for research purposes, please cite the above paper
// in any resulting publication.
//
// EXAMPLE
//
// Generate point samples whose density is guided by a 512x512 grayscale image:
//
//     int npoints;
//     float* points;
//     int maxpoints = 1e6;
//     float density = 30000;
//     par_bluenoise_context* ctx;
//     ctx = par_bluenoise_from_file("bluenoise.bin", maxpoints);
//     par_bluenoise_density_from_gray(ctx, source_pixels, 512, 512, 1);
//     points = par_bluenoise_generate(ctx, density, &npoints);
//     ... Draw points here.  Each point is a three-tuple of (X Y RANK).
//     par_bluenoise_free(ctx);
//
// Distributed under the MIT License, see bottom of file.

#ifndef PAR_BLUENOISE_H
#define PAR_BLUENOISE_H

#ifdef __cplusplus
extern "C" {
#endif

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

typedef unsigned char par_byte;

// Encapsulates a tile set and an optional density function.
typedef struct par_bluenoise_context_s par_bluenoise_context;

// Creates a bluenoise context using the given tileset.  The first argument is
// the file path the bin file.  The second argument is the maximum number of
// points that you expect to ever be generated.
par_bluenoise_context* par_bluenoise_from_file(char const* path, int maxpts);

// Creates a bluenoise context using the given tileset.  The first and second
// arguments describe a memory buffer containing the contents of the bin file.
// The third argument is the maximum number of points that you expect to ever
// be generated.
par_bluenoise_context* par_bluenoise_from_buffer(
    par_byte const* buffer, int nbytes, int maxpts);

// Sets up a scissoring rectangle using the given lower-left and upper-right
// coordinates.  By default the scissor encompasses [-0.5, -0.5] - [0.5, 0.5],
// which is the entire sampling domain for the two "generate" methods.
void par_bluenoise_set_viewport(
    par_bluenoise_context*, float left, float bottom, float right, float top);

// Sets up a reference window size.  The only purpose of this is to ensure
// that apparent density remains constant when the window gets resized.
// Clients should call this *before* calling par_bluenoise_set_viewport.
void par_bluenoise_set_window(par_bluenoise_context*, int width, int height);

// Frees all memory associated with the given bluenoise context.
void par_bluenoise_free(par_bluenoise_context* ctx);

// Copies a grayscale image into the bluenoise context to guide point density.
// Darker regions generate a higher number of points. The given bytes-per-pixel
// value is the stride between pixels.
void par_bluenoise_density_from_gray(par_bluenoise_context* ctx,
    const unsigned char* pixels, int width, int height, int bpp);

// Creates a binary mask to guide point density. The given bytes-per-pixel
// value is the stride between pixels, which must be 4 or less.
void par_bluenoise_density_from_color(par_bluenoise_context* ctx,
    const unsigned char* pixels, int width, int height, int bpp,
    unsigned int background_color, int invert);

// Generates samples using Recursive Wang Tiles.  This is really fast!
// The returned pointer is a list of three-tuples, where XY are in [-0.5, +0.5]
// and Z is a rank value that can be used to create a progressive ordering.
// The caller should not free the returned pointer.
float* par_bluenoise_generate(
    par_bluenoise_context* ctx, float density, int* npts);

// Generates an ordered sequence of tuples with the specified sequence length.
// This is slower than the other "generate" method because it uses a dumb
// backtracking method to determine a reasonable density value, and it
// automatically sorts the output by rank.  The dims argument must be 2 or more;
// it represents the desired stride (in floats) between consecutive verts in the
// returned data buffer.
float* par_bluenoise_generate_exact(
    par_bluenoise_context* ctx, int npts, int dims);

// Performs an in-place sort of 3-tuples, based on the 3rd component, then
// replaces the 3rd component with an index.
void par_bluenoise_sort_by_rank(float* pts, int npts);

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

#ifdef PAR_BLUENOISE_IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#define PAR_MINI(a, b) ((a < b) ? a : b)
#define PAR_MAXI(a, b) ((a > b) ? a : b)

typedef struct {
    float x;
    float y;
} par_vec2;

typedef struct {
    float x;
    float y;
    float rank;
} par_vec3;

typedef struct {
    int n, e, s, w;
    int nsubtiles, nsubdivs, npoints, nsubpts;
    int** subdivs;
    par_vec2* points;
    par_vec2* subpts;
} par_tile;

struct par_bluenoise_context_s {
    par_vec3* points;
    par_tile* tiles;
    float global_density;
    float left, bottom, right, top;
    int ntiles, nsubtiles, nsubdivs;
    int npoints;
    int maxpoints;
    int density_width;
    int density_height;
    unsigned char* density;
    float mag;
    int window_width;
    int window_height;
    int abridged;
};

static float sample_density(par_bluenoise_context* ctx, float x, float y)
{
    unsigned char* density = ctx->density;
    if (!density) {
        return 1;
    }
    int width = ctx->density_width;
    int height = ctx->density_height;
    y = 1 - y;
    x -= 0.5;
    y -= 0.5;
    float tx = x * PAR_MAXI(width, height);
    float ty = y * PAR_MAXI(width, height);
    x += 0.5;
    y += 0.5;
    tx += width / 2;
    ty += height / 2;
    int ix = PAR_CLAMP((int) tx, 0, width - 2);
    int iy = PAR_CLAMP((int) ty, 0, height - 2);
    return density[iy * width + ix] / 255.0f;
}

static void recurse_tile(
    par_bluenoise_context* ctx, par_tile* tile, float x, float y, int level)
{
    float left = ctx->left, right = ctx->right;
    float top = ctx->top, bottom = ctx->bottom;
    float mag = ctx->mag;
    float tileSize = 1.f / powf(ctx->nsubtiles, level);
    if (x + tileSize < left || x > right || y + tileSize < bottom || y > top) {
        return;
    }
    float depth = powf(ctx->nsubtiles, 2 * level);
    float threshold = mag / depth * ctx->global_density - tile->npoints;
    int ntests = PAR_MINI(tile->nsubpts, threshold);
    float factor = 1.f / mag * depth / ctx->global_density;
    for (int i = 0; i < ntests; i++) {
        float px = x + tile->subpts[i].x * tileSize;
        float py = y + tile->subpts[i].y * tileSize;
        if (px < left || px > right || py < bottom || py > top) {
            continue;
        }
        if (sample_density(ctx, px, py) < (i + tile->npoints) * factor) {
            continue;
        }
        ctx->points[ctx->npoints].x = px - 0.5;
        ctx->points[ctx->npoints].y = py - 0.5;
        ctx->points[ctx->npoints].rank = (level + 1) + i * factor;
        ctx->npoints++;
        if (ctx->npoints >= ctx->maxpoints) {
            return;
        }
    }
    const float scale = tileSize / ctx->nsubtiles;
    if (threshold <= tile->nsubpts) {
        return;
    }
    level++;
    for (int ty = 0; ty < ctx->nsubtiles; ty++) {
        for (int tx = 0; tx < ctx->nsubtiles; tx++) {
            int tileIndex = tile->subdivs[0][ty * ctx->nsubtiles + tx];
            par_tile* subtile = &ctx->tiles[tileIndex];
            recurse_tile(ctx, subtile, x + tx * scale, y + ty * scale, level);
        }
    }
}

void par_bluenoise_set_window(par_bluenoise_context* ctx, int width, int height)
{
    ctx->window_width = width;
    ctx->window_height = height;
}

void par_bluenoise_set_viewport(par_bluenoise_context* ctx, float left,
    float bottom, float right, float top)
{
    // Transform [-.5, +.5] to [0, 1]
    left = ctx->left = left + 0.5;
    right = ctx->right = right + 0.5;
    bottom = ctx->bottom = bottom + 0.5;
    top = ctx->top = top + 0.5;

    // Determine magnification factor BEFORE clamping.
    float scale = 1000 * (top - bottom) / ctx->window_height;
    ctx->mag = powf(scale, -2);

    // The density function is only sampled in [0, +1].
    ctx->left = PAR_CLAMP(left, 0, 1);
    ctx->right = PAR_CLAMP(right, 0, 1);
    ctx->bottom = PAR_CLAMP(bottom, 0, 1);
    ctx->top = PAR_CLAMP(top, 0, 1);
}

float* par_bluenoise_generate(
    par_bluenoise_context* ctx, float density, int* npts)
{
    ctx->global_density = density;
    ctx->npoints = 0;
    float left = ctx->left;
    float right = ctx->right;
    float bottom = ctx->bottom;
    float top = ctx->top;
    float mag = ctx->mag;

    int ntests = PAR_MINI(ctx->tiles[0].npoints, mag * ctx->global_density);
    float factor = 1.f / mag / ctx->global_density;
    for (int i = 0; i < ntests; i++) {
        float px = ctx->tiles[0].points[i].x;
        float py = ctx->tiles[0].points[i].y;
        if (px < left || px > right || py < bottom || py > top) {
            continue;
        }
        if (sample_density(ctx, px, py) < (i + 1) * factor) {
            continue;
        }
        ctx->points[ctx->npoints].x = px - 0.5;
        ctx->points[ctx->npoints].y = py - 0.5;
        ctx->points[ctx->npoints].rank = i * factor;
        ctx->npoints++;
        if (ctx->npoints >= ctx->maxpoints) {
            break;
        }
    }
    recurse_tile(ctx, &ctx->tiles[0], 0, 0, 0);
    *npts = ctx->npoints;
    return &ctx->points->x;
}

#define freadi()   \
    *((int*) ptr); \
    ptr += sizeof(int)

#define freadf()     \
    *((float*) ptr); \
    ptr += sizeof(float)

static par_bluenoise_context* par_bluenoise_create(
    char const* filepath, int nbytes, int maxpts)
{
    par_bluenoise_context* ctx = PAR_MALLOC(par_bluenoise_context, 1);
    ctx->maxpoints = maxpts;
    ctx->points = PAR_MALLOC(par_vec3, maxpts);
    ctx->density = 0;
    ctx->abridged = 0;
    par_bluenoise_set_window(ctx, 1024, 768);
    par_bluenoise_set_viewport(ctx, -.5, -.5, .5, .5);

    char* buf = 0;
    if (nbytes == 0) {
        FILE* fin = fopen(filepath, "rb");
        assert(fin);
        fseek(fin, 0, SEEK_END);
        nbytes = (int) ftell(fin);
        fseek(fin, 0, SEEK_SET);
        buf = PAR_MALLOC(char, nbytes);
        int consumed = (int) fread(buf, nbytes, 1, fin);
        assert(consumed == 1);
        fclose(fin);
    }

    char const* ptr = buf ? buf : filepath;
    int ntiles = ctx->ntiles = freadi();
    int nsubtiles = ctx->nsubtiles = freadi();
    int nsubdivs = ctx->nsubdivs = freadi();
    par_tile* tiles = ctx->tiles = PAR_MALLOC(par_tile, ntiles);
    for (int i = 0; i < ntiles; i++) {
        tiles[i].n = freadi();
        tiles[i].e = freadi();
        tiles[i].s = freadi();
        tiles[i].w = freadi();
        tiles[i].subdivs = PAR_MALLOC(int*, nsubdivs);
        for (int j = 0; j < nsubdivs; j++) {
            int* subdiv = PAR_MALLOC(int, PAR_SQR(nsubtiles));
            for (int k = 0; k < PAR_SQR(nsubtiles); k++) {
                subdiv[k] = freadi();
            }
            tiles[i].subdivs[j] = subdiv;
        }
        tiles[i].npoints = freadi();
        tiles[i].points = PAR_MALLOC(par_vec2, tiles[i].npoints);
        for (int j = 0; j < tiles[i].npoints; j++) {
            tiles[i].points[j].x = freadf();
            tiles[i].points[j].y = freadf();
        }
        tiles[i].nsubpts = freadi();
        tiles[i].subpts = PAR_MALLOC(par_vec2, tiles[i].nsubpts);
        for (int j = 0; j < tiles[i].nsubpts; j++) {
            tiles[i].subpts[j].x = freadf();
            tiles[i].subpts[j].y = freadf();
        }

        // The following hack allows for an optimization whereby
        // the first tile's point set is re-used for every other tile.
        // This goes against the entire purpose of Recursive Wang Tiles,
        // but in many applications the qualatitive loss is not
        // observable, and the footprint savings are huge (10x).

        if (tiles[i].npoints == 0) {
            ctx->abridged = 1;
            tiles[i].npoints = tiles[0].npoints;
            tiles[i].points = tiles[0].points;
            tiles[i].nsubpts = tiles[0].nsubpts;
            tiles[i].subpts = tiles[0].subpts;
        }
    }
    free(buf);
    return ctx;
}

par_bluenoise_context* par_bluenoise_from_file(char const* path, int maxpts)
{
    return par_bluenoise_create(path, 0, maxpts);
}

par_bluenoise_context* par_bluenoise_from_buffer(
    par_byte const* buffer, int nbytes, int maxpts)
{
    return par_bluenoise_create((char const*) buffer, nbytes, maxpts);
}

void par_bluenoise_density_from_gray(par_bluenoise_context* ctx,
    const unsigned char* pixels, int width, int height, int bpp)
{
    ctx->density_width = width;
    ctx->density_height = height;
    ctx->density = PAR_MALLOC(unsigned char, width * height);
    unsigned char* dst = ctx->density;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            *dst++ = 255 - (*pixels);
            pixels += bpp;
        }
    }
}

void par_bluenoise_density_from_color(par_bluenoise_context* ctx,
    const unsigned char* pixels, int width, int height, int bpp,
    unsigned int background_color, int invert)
{
    unsigned int bkgd = background_color;
    ctx->density_width = width;
    ctx->density_height = height;
    ctx->density = PAR_MALLOC(unsigned char, width * height);
    unsigned char* dst = ctx->density;
    unsigned int mask = 0x000000ffu;
    if (bpp > 1) {
        mask |= 0x0000ff00u;
    }
    if (bpp > 2) {
        mask |= 0x00ff0000u;
    }
    if (bpp > 3) {
        mask |= 0xff000000u;
    }
    assert(bpp <= 4);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            unsigned int val = (*((unsigned int*) pixels)) & mask;
            val = invert ? (val == bkgd) : (val != bkgd);
            *dst++ = val ? 255 : 0;
            pixels += bpp;
        }
    }
}

void par_bluenoise_free(par_bluenoise_context* ctx)
{
    free(ctx->points);
    for (int t = 0; t < ctx->ntiles; t++) {
        for (int s = 0; s < ctx->nsubdivs; s++) {
            free(ctx->tiles[t].subdivs[s]);
        }
        free(ctx->tiles[t].subdivs);
        if (t == 0 || !ctx->abridged) {
            free(ctx->tiles[t].points);
            free(ctx->tiles[t].subpts);
        }
    }
    free(ctx->tiles);
    free(ctx->density);
}

int cmp(const void* a, const void* b)
{
    const par_vec3* v1 = (const par_vec3*) a;
    const par_vec3* v2 = (const par_vec3*) b;
    if (v1->rank < v2->rank) {
        return -1;
    }
    if (v1->rank > v2->rank) {
        return 1;
    }
    return 0;
}

void par_bluenoise_sort_by_rank(float* floats, int npts)
{
    par_vec3* vecs = (par_vec3*) floats;
    qsort(vecs, npts, sizeof(vecs[0]), cmp);
    for (int i = 0; i < npts; i++) {
        vecs[i].rank = i;
    }
}

float* par_bluenoise_generate_exact(
    par_bluenoise_context* ctx, int npts, int stride)
{
    assert(stride >= 2);
    int maxpoints = npts * 2;
    if (ctx->maxpoints < maxpoints) {
        free(ctx->points);
        ctx->maxpoints = maxpoints;
        ctx->points = PAR_MALLOC(par_vec3, maxpoints);
    }
    int ngenerated = 0;
    int nprevious = 0;
    int ndesired = npts;
    float density = 2048;
    while (ngenerated < ndesired) {
        par_bluenoise_generate(ctx, density, &ngenerated);

        // Might be paranoid, but break if something fishy is going on:
        if (ngenerated == nprevious) {
            return 0;
        }

        // Perform crazy heuristic to approach a nice density:
        if (ndesired / ngenerated >= 2) {
            density *= 2;
        } else {
            density += density / 10;
        }

        nprevious = ngenerated;
    }
    par_bluenoise_sort_by_rank(&ctx->points->x, ngenerated);
    if (stride != 3) {
        int nbytes = sizeof(float) * stride * ndesired;
        float* pts = PAR_MALLOC(float, stride * ndesired);
        float* dst = pts;
        const float* src = &ctx->points->x;
        for (int i = 0; i < ndesired; i++, src++) {
            *dst++ = *src++;
            *dst++ = *src++;
            if (stride > 3) {
                *dst++ = *src;
                dst += stride - 3;
            }
        }
        memcpy(ctx->points, pts, nbytes);
        free(pts);
    }
    return &ctx->points->x;
}

#undef PAR_MINI
#undef PAR_MAXI

#endif // PAR_BLUENOISE_IMPLEMENTATION
#endif // PAR_BLUENOISE_H

// par_bluenoise is distributed under the MIT license:
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
