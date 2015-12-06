// MSQUARES :: https://github.com/prideout/par
// Converts fp32 grayscale images, or 8-bit color images, into triangles.
//
// For grayscale images, a threshold is specified to determine insideness.
// Color images can be r8, rg16, rgb24, or rgba32.  In all cases, the exact
// color of each pixel determines insideness.  For a visual overview of
// the API, see:
//
//     https://....
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

typedef uint8_t par_byte;

typedef struct par_msquares_meshlist_s par_msquares_meshlist;

// Encapsulates the results of a marching squares operation.
typedef struct {
    float* points;        // pointer to XY (or XYZ) vertex coordinates
    int npoints;          // number of vertex coordinates
    uint16_t* triangles;  // pointer to 3-tuples of vertex indices
    int ntriangles;       // number of 3-tuples
    int dim;              // number of floats per point (either 2 or 3)
} par_msquares_mesh;

// Reverses the "insideness" test.
// Only applicable with from_grayscale and from_color.
#define PAR_MSQUARES_INVERT (1 << 0)

// Returns a meshlist with two meshes: one for the inside, one for the outside.
// Only applicable with from_grayscale and from_color.
#define PAR_MSQUARES_DUAL (1 << 1)

// Returned meshes have 3-tuple coordinates instead of 2-tuples.
// With from_color and from_colors, the Z coordinate represents the alpha value
// of the color.  With from_grayscale and from_levels, the Z coordinate
// represents the value of the nearest pixel in the source image.
#define PAR_MSQUARES_HEIGHTS (1 << 2)

// Adds extrusion triangles to each mesh other than the lowest mesh.  Requires
// the PAR_MSQUARES_HEIGHTS flag to be present.
#define PAR_MSQUARES_CONNECT (1 << 3)

// Enables quicky & dirty (not maximal) simpification of the returned mesh.
#define PAR_MSQUARES_SIMPLIFY (1 << 4)

par_msquares_meshlist* par_msquares_from_grayscale(float const* data, int width,
    int height, int cellsize, float threshold, int flags);

par_msquares_meshlist* par_msquares_from_levels(float const* data, int width,
    int height, int cellsize, float const* thresholds, int nthresholds,
    int flags);

par_msquares_meshlist* par_msquares_from_color(par_byte const* data, int width,
    int height, int cellsize, uint32_t color, int bpp, int flags);

par_msquares_meshlist* par_msquares_from_colors(par_byte const* data, int width,
    int height, int cellsize, uint32_t const* colors, int ncolors, int bpp,
    int flags);

typedef int (*par_msquares_inside_fn)(int, void*);
typedef float (*par_msquares_height_fn)(float, float, void*);

par_msquares_meshlist* par_msquares_from_function(int width, int height,
    int cellsize, int flags, void* context, par_msquares_inside_fn insidefn,
    par_msquares_height_fn heightfn);

par_msquares_mesh* par_msquares_get_mesh(par_msquares_meshlist*, int n);

int par_msquares_get_count(par_msquares_meshlist*);

void par_msquares_free(par_msquares_meshlist*);

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_MSQUARES_IMPLEMENTATION

#define MIN(a, b) (a > b ? b : a)
#define MAX(a, b) (a > b ? a : b)
#define CLAMP(v, lo, hi) MAX(lo, MIN(hi, v))

struct par_msquares_meshlist_s {
    int nmeshes;
    par_msquares_mesh** meshes;
};

static int** point_table = 0;
static int** triangle_table = 0;

static void par_init_tables()
{
    point_table = malloc(16 * sizeof(int*));
    triangle_table = malloc(16 * sizeof(int*));

    char const* CODE_TABLE =
        "0 0\n"
        "1 1 0 1 7\n"
        "2 1 1 2 3\n"
        "3 2 0 2 3 3 7 0\n"
        "4 1 7 5 6\n"
        "5 2 0 1 5 5 6 0\n"
        "6 2 1 2 3 7 5 6\n"
        "7 3 0 2 3 0 3 5 0 5 6\n"
        "8 1 3 4 5\n"
        "9 4 0 1 3 0 3 4 0 4 5 0 5 7\n"
        "a 2 1 2 4 4 5 1\n"
        "b 3 0 2 4 0 4 5 0 5 7\n"
        "c 2 7 3 4 4 6 7\n"
        "d 3 0 1 3 0 3 4 0 4 6\n"
        "e 3 1 2 4 1 4 6 1 6 7\n"
        "f 2 0 2 4 4 6 0\n";

    char const* table_token = CODE_TABLE;
    for (int i = 0; i < 16; i++) {
        char code = *table_token;
        assert(i == code - (code >= 'a' ? ('a' - 0xa) : '0'));
        table_token += 2;
        int ntris = *table_token - '0';
        table_token += 2;
        int* sqrtris = triangle_table[i] =
                malloc((ntris + 1) * 3 * sizeof(int));
        sqrtris[0] = ntris;
        int mask = 0;
        int* sqrpts = point_table[i] = malloc(7 * sizeof(int));
        sqrpts[0] = 0;
        for (int j = 0; j < ntris * 3; j++, table_token += 2) {
            int midp = *table_token - '0';
            int bit = 1 << midp;
            if (!(mask & bit)) {
                mask |= bit;
                sqrpts[++sqrpts[0]] = midp;
            }
            sqrtris[j + 1] = midp;
        }
    }
}

typedef struct {
    float const* data;
    float threshold;
    int width;
    int height;
} gray_context;

static int gray_inside(int location, void* contextptr)
{
    gray_context* context = (gray_context*) contextptr;
    return context->data[location] > context->threshold;
}

static float gray_height(float x, float y, void* contextptr)
{
    gray_context* context = (gray_context*) contextptr;
    int i = CLAMP(context->width * x, 0, context->width - 1);
    int j = CLAMP(context->height * y, 0, context->height - 1);
    return context->data[i + j * context->width];
}

typedef struct {
    par_byte const* data;
    par_byte color[4];
    int bpp;
    int width;
    int height;
} color_context;

static int color_inside(int location, void* contextptr)
{
    color_context* context = (color_context*) contextptr;
    par_byte const* data = context->data + location * context->bpp;
    for (int i = 0; i < context->bpp; i++) {
        if (data[i] != context->color[i]) {
            return 0;
        }
    }
    return 1;
}

static float color_height(float x, float y, void* contextptr)
{
    color_context* context = (color_context*) contextptr;
    assert(context->bpp == 4);
    int i = CLAMP(context->width * x, 0, context->width - 1);
    int j = CLAMP(context->height * y, 0, context->height - 1);
    int k = i + j * context->width;
    return context->data[k * 4 + 3] / 255.0;
}

par_msquares_meshlist* par_msquares_from_color(par_byte const* data, int width,
    int height, int cellsize, uint32_t color, int bpp, int flags)
{
    color_context context;
    context.bpp = bpp;
    context.color[0] = (color >> 16) & 0xff;
    context.color[1] = (color >> 8) & 0xff;
    context.color[2] = (color & 0xff);
    context.color[3] = (color >> 24) & 0xff;
    context.data = data;
    context.width = width;
    context.height = height;
    return par_msquares_from_function(
        width, height, cellsize, flags, &context, color_inside, color_height);
}

par_msquares_meshlist* par_msquares_from_grayscale(float const* data, int width,
    int height, int cellsize, float threshold, int flags)
{
    gray_context context;
    context.width = width;
    context.height = height;
    context.data = data;
    context.threshold = threshold;
    return par_msquares_from_function(
        width, height, cellsize, flags, &context, gray_inside, gray_height);
}

par_msquares_mesh* par_msquares_get_mesh(
    par_msquares_meshlist* mlist, int mindex)
{
    assert(mlist && mindex < mlist->nmeshes);
    return mlist->meshes[mindex];
}

int par_msquares_get_count(par_msquares_meshlist* mlist)
{
    assert(mlist);
    return mlist->nmeshes;
}

void par_msquares_free(par_msquares_meshlist* mlist)
{
    par_msquares_mesh** meshes = mlist->meshes;
    for (int i = 0; i < mlist->nmeshes; i++) {
        free(meshes[i]->points);
        free(meshes[i]->triangles);
        free(meshes[i]);
    }
    free(meshes);
    free(mlist);
}

static par_msquares_meshlist* par_msquares_merge(par_msquares_meshlist** lists,
    int count)
{
    par_msquares_meshlist* merged = malloc(sizeof(par_msquares_meshlist));
    merged->nmeshes = 0;
    for (int i = 0; i < count; i++) {
        merged->nmeshes += lists[i]->nmeshes;
    }
    merged->meshes = malloc(sizeof(par_msquares_mesh*) * merged->nmeshes);
    par_msquares_mesh** pmesh = merged->meshes;
    for (int i = 0; i < count; i++) {
        par_msquares_meshlist* meshlist = lists[i];
        for (int j = 0; j < meshlist->nmeshes; j++) {
            *pmesh++ = meshlist->meshes[j];
        }
        free(meshlist);
    }
    return merged;
}

par_msquares_meshlist* par_msquares_from_function(int width, int height,
    int cellsize, int flags, void* context, par_msquares_inside_fn insidefn,
    par_msquares_height_fn heightfn)
{
    assert(width > 0 && width % cellsize == 0);
    assert(height > 0 && height % cellsize == 0);

    if (flags & PAR_MSQUARES_DUAL) {
        flags &= ~PAR_MSQUARES_INVERT;
        flags &= ~PAR_MSQUARES_DUAL;
        par_msquares_meshlist* m[2];
        m[0] = par_msquares_from_function(width, height, cellsize, flags,
            context, insidefn, heightfn);
        flags |= PAR_MSQUARES_INVERT;
        m[1] = par_msquares_from_function(width, height, cellsize, flags,
            context, insidefn, heightfn);
        return par_msquares_merge(m, 2);
    }

    int invert = flags & PAR_MSQUARES_INVERT;

    // Create the two code tables if we haven't already.  These are tables of
    // fixed constants, so it's embarassing that we use dynamic memory
    // allocation for them.  However it's easy and it's one-time-only.

    if (!point_table) {
        par_init_tables();
    }

    // Allocate the meshlist and the first mesh.

    par_msquares_meshlist* mlist = malloc(sizeof(par_msquares_meshlist));
    mlist->nmeshes = 1;
    mlist->meshes = malloc(sizeof(par_msquares_mesh*));
    mlist->meshes[0] = malloc(sizeof(par_msquares_mesh));
    par_msquares_mesh* mesh = mlist->meshes[0];
    mesh->dim = (flags & PAR_MSQUARES_HEIGHTS) ? 3 : 2;
    int ncols = width / cellsize;
    int nrows = height / cellsize;

    // Worst case is four triangles and six verts per cell, so allocate that
    // much.

    int maxtris = ncols * nrows * 4;
    uint16_t* tris = malloc(maxtris * 3 * sizeof(uint16_t));
    int ntris = 0;
    int maxpts = ncols * nrows * 6;
    float* pts = malloc(maxpts * mesh->dim * sizeof(float));
    int npts = 0;

    // The "verts" x/y/z arrays are the 4 corners and 4 midpoints around the
    // square,
    // in counter-clockwise order.  The origin of "triangle space" is at the
    // lower-left, although we expect the image data to be in raster order
    // (starts at top-left).

    float normalization = 1.0f / MAX(width, height);
    float normalized_cellsize = cellsize * normalization;
    int maxrow = (height - 1) * width;
    uint16_t* ptris = tris;
    float* ppts = pts;
    float vertsx[8], vertsy[8];
    int* prevrowmasks = calloc(sizeof(int) * ncols, 1);
    int* prevrowinds = calloc(sizeof(int) * ncols * 3, 1);

    // Do the march!

    for (int row = 0; row < nrows; row++) {
        vertsx[0] = vertsx[6] = vertsx[7] = 0;
        vertsx[1] = vertsx[5] = 0.5 * normalized_cellsize;
        vertsx[2] = vertsx[3] = vertsx[4] = normalized_cellsize;
        vertsy[0] = vertsy[1] = vertsy[2] = normalized_cellsize * (row + 1);
        vertsy[4] = vertsy[5] = vertsy[6] = normalized_cellsize * row;
        vertsy[3] = vertsy[7] = normalized_cellsize * (row + 0.5);

        int northi = row * cellsize * width;
        int southi = MIN(northi + cellsize * width, maxrow);
        int northwest = invert ^ insidefn(northi, context);
        int southwest = invert ^ insidefn(southi, context);
        int previnds[8] = {0};
        int prevmask = 0;

        for (int col = 0; col < ncols; col++) {
            northi += cellsize;
            southi += cellsize;
            if (col == ncols - 1) {
                northi--;
                southi--;
            }

            int northeast = invert ^ insidefn(northi, context);
            int southeast = invert ^ insidefn(southi, context);
            int code = southwest | (southeast << 1) | (northwest << 2) |
                (northeast << 3);

            int const* pointspec = point_table[code];
            int ptspeclength = *pointspec++;
            int currinds[8] = {0};
            int mask = 0;
            int prevrowmask = prevrowmasks[col];
            while (ptspeclength--) {
                int midp = *pointspec++;
                int bit = 1 << midp;
                mask |= bit;

                // The following six conditionals perform welding to reduce the
                // number of vertices.  The first three perform welding with the
                // cell to the west; the latter three perform welding with the
                // cell to the north.

                if (bit == 1 && (prevmask & 4)) {
                    currinds[midp] = previnds[2];
                    continue;
                }
                if (bit == 128 && (prevmask & 8)) {
                    currinds[midp] = previnds[3];
                    continue;
                }
                if (bit == 64 && (prevmask & 16)) {
                    currinds[midp] = previnds[4];
                    continue;
                }
                if (bit == 16 && (prevrowmask & 4)) {
                    currinds[midp] = prevrowinds[col * 3 + 2];
                    continue;
                }
                if (bit == 32 && (prevrowmask & 2)) {
                    currinds[midp] = prevrowinds[col * 3 + 1];
                    continue;
                }
                if (bit == 64 && (prevrowmask & 1)) {
                    currinds[midp] = prevrowinds[col * 3 + 0];
                    continue;
                }

                ppts[0] = vertsx[midp];
                ppts[1] = vertsy[midp];

                // Adjust the midpoints to a more exact crossing point.
                if (midp == 1) {
                    int begin = southi - cellsize / 2;
                    int previous = 0;
                    for (int i = 0; i < cellsize; i++) {
                        int offset = begin + i / 2 * ((i % 2) ? -1 : 1);
                        int inside = insidefn(offset, context);
                        if (i > 0 && inside != previous) {
                            ppts[0] = normalization *
                                (col * cellsize + offset - southi + cellsize);
                            break;
                        }
                        previous = inside;
                    }
                } else if (midp == 5) {
                    int begin = northi - cellsize / 2;
                    int previous = 0;
                    for (int i = 0; i < cellsize; i++) {
                        int offset = begin + i / 2 * ((i % 2) ? -1 : 1);
                        int inside = insidefn(offset, context);
                        if (i > 0 && inside != previous) {
                            ppts[0] = normalization *
                                (col * cellsize + offset - northi + cellsize);
                            break;
                        }
                        previous = inside;
                    }
                } else if (midp == 3) {
                    int begin = northi + width * cellsize / 2;
                    int previous = 0;
                    for (int i = 0; i < cellsize; i++) {
                        int offset = begin +
                            width * (i / 2 * ((i % 2) ? -1 : 1));
                        int inside = insidefn(offset, context);
                        if (i > 0 && inside != previous) {
                            ppts[1] = normalization *
                                (row * cellsize +
                                (offset - northi) / (float) width);
                            break;
                        }
                        previous = inside;
                    }
                } else if (midp == 7) {
                    int begin = northi + width * cellsize / 2 - cellsize;
                    int previous = 0;
                    for (int i = 0; i < cellsize; i++) {
                        int offset = begin +
                            width * (i / 2 * ((i % 2) ? -1 : 1));
                        int inside = insidefn(offset, context);
                        if (i > 0 && inside != previous) {
                            ppts[1] = normalization *
                                (row * cellsize +
                                (offset - northi - cellsize) / (float) width);
                            break;
                        }
                        previous = inside;
                    }
                }

                if (mesh->dim == 3) {
                    ppts[2] = heightfn(ppts[0], ppts[1], context);
                }

                ppts += mesh->dim;
                currinds[midp] = npts++;
            }

            // Add triangles.
            int const* trianglespec = triangle_table[code];
            int trispeclength = *trianglespec++;
            while (trispeclength--) {
                int a = *trianglespec++;
                int b = *trianglespec++;
                int c = *trianglespec++;
                *ptris++ = currinds[c];
                *ptris++ = currinds[b];
                *ptris++ = currinds[a];
                ntris++;
            }

            // Prepare for the next cell.
            prevrowmasks[col] = mask;
            prevrowinds[col * 3 + 0] = currinds[0];
            prevrowinds[col * 3 + 1] = currinds[1];
            prevrowinds[col * 3 + 2] = currinds[2];
            prevmask = mask;
            northwest = northeast;
            southwest = southeast;
            for (int i = 0; i < 8; i++) {
                previnds[i] = currinds[i];
                vertsx[i] += normalized_cellsize;
            }
        }
    }

    free(prevrowmasks);
    free(prevrowinds);
    assert(npts <= maxpts);
    assert(ntris <= maxtris);
    mesh->npoints = npts;
    mesh->points = pts;
    mesh->ntriangles = ntris;
    mesh->triangles = tris;
    return mlist;
}

#undef MIN
#undef MAX
#undef CLAMP
#endif
