// MSQUARES :: https://github.com/prideout/par
// Converts fp32 grayscale images, or 8-bit color images, into triangles.
//
// For grayscale images, a threshold is specified to determine insideness.
// For color images, an exact color is specified to determine insideness.
// Color images can be r8, rg16, rgb24, or rgba32. For a visual overview of
// the API and all the flags, see:
//
//     http://github.prideout.net/marching-squares/
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

#include <stdint.h>

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
    int nconntriangles;   // internal use only
} par_msquares_mesh;

// Reverses the "insideness" test.
#define PAR_MSQUARES_INVERT (1 << 0)

// Returns a meshlist with two meshes: one for the inside, one for the outside.
#define PAR_MSQUARES_DUAL (1 << 1)

// Returned meshes have 3-tuple coordinates instead of 2-tuples. When using
// from_color, the Z coordinate represents the alpha value of the color.  With
// from_grayscale, the Z coordinate represents the value of the nearest pixel in
// the source image.
#define PAR_MSQUARES_HEIGHTS (1 << 2)

// Applies a step function to the Z coordinates.  Requires HEIGHTS and DUAL.
#define PAR_MSQUARES_SNAP (1 << 3)

// Adds extrusion triangles to each mesh other than the lowest mesh.  Requires
// the PAR_MSQUARES_HEIGHTS flag to be present.
#define PAR_MSQUARES_CONNECT (1 << 4)

// Enables quick & dirty (not best) simpification of the returned mesh.
#define PAR_MSQUARES_SIMPLIFY (1 << 5)

par_msquares_meshlist* par_msquares_grayscale(float const* data, int width,
    int height, int cellsize, float threshold, int flags);

par_msquares_meshlist* par_msquares_color(par_byte const* data, int width,
    int height, int cellsize, uint32_t color, int bpp, int flags);

par_msquares_mesh const* par_msquares_get_mesh(par_msquares_meshlist*, int n);

int par_msquares_get_count(par_msquares_meshlist*);

void par_msquares_free(par_msquares_meshlist*);

typedef int (*par_msquares_inside_fn)(int, void*);
typedef float (*par_msquares_height_fn)(float, float, void*);

par_msquares_meshlist* par_msquares_function(int width, int height,
    int cellsize, int flags, void* context, par_msquares_inside_fn insidefn,
    par_msquares_height_fn heightfn);

par_msquares_meshlist* par_msquares_grayscale_multi(float const* data,
    int width, int height, int cellsize, float const* thresholds,
    int nthresholds, int flags);

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_MSQUARES_IMPLEMENTATION

#include <stdlib.h>
#include <assert.h>
#include <float.h>

#define PAR_MIN(a, b) (a > b ? b : a)
#define PAR_MAX(a, b) (a > b ? a : b)
#define PAR_CLAMP(v, lo, hi) PAR_MAX(lo, PAR_MIN(hi, v))
#define PAR_ALLOC(T, N) ((T*) calloc(N * sizeof(T), 1))

struct par_msquares_meshlist_s {
    int nmeshes;
    par_msquares_mesh** meshes;
};

static int** par_msquares_binary_point_table = 0;
static int** par_msquares_binary_triangle_table = 0;
static int** par_msquares_quaternary_triangle_table = 0;

static par_msquares_meshlist* par_msquares_merge(par_msquares_meshlist** lists,
    int count, int snap);

static void par_init_tables()
{
    char const* BINARY_TABLE =
        "0"
        "1017"
        "1123"
        "2023370"
        "1756"
        "2015560"
        "2123756"
        "3023035056"
        "1345"
        "4013034045057"
        "2124451"
        "3024045057"
        "2734467"
        "3013034046"
        "3124146167"
        "2024460";
    char const* binary_token = BINARY_TABLE;

    par_msquares_binary_point_table = PAR_ALLOC(int*, 16);
    par_msquares_binary_triangle_table = PAR_ALLOC(int*, 16);
    for (int i = 0; i < 16; i++) {
        int ntris = *binary_token - '0';
        binary_token++;
        int* sqrtris = par_msquares_binary_triangle_table[i] =
            PAR_ALLOC(int, (ntris + 1) * 3);
        sqrtris[0] = ntris;
        int mask = 0;
        int* sqrpts = par_msquares_binary_point_table[i] = PAR_ALLOC(int, 7);
        sqrpts[0] = 0;
        for (int j = 0; j < ntris * 3; j++, binary_token++) {
            int midp = *binary_token - '0';
            int bit = 1 << midp;
            if (!(mask & bit)) {
                mask |= bit;
                sqrpts[++sqrpts[0]] = midp;
            }
            sqrtris[j + 1] = midp;
        }
    }

    char const* QUATERNARY_TABLE =
        "2024046000"
        "3702724745175600"
        "3702724745017560"
        "3702724745001756"
        "3560502523134500"
        "2023037273474600"
        "2023037283484527857560"
        "2023037283484502785756"
        "3560502523013450"
        "2023037278575628348450"
        "2023037027347460"
        "2023037028348452785756"
        "3560502523001345"
        "2023037278575602834845"
        "2023037027857562834845"
        "2023037002734746"
        "3346360301112300"
        "4018087834845412313878575600"
        "4013034045057112317560"
        "4013034045057112301756"
        "2015056212414500"
        "1701312414616700"
        "2018087312414616727857560"
        "2018087312414616702785756"
        "2015056212313828348450"
        "1701467161262363513450"
        "2018087212313827347460"
        "2018087212313828348452785756"
        "2015056212313802834845"
        "1701467161262363501345"
        "2018087212313827857562834845"
        "2018087212313802734746"
        "3346360301011230"
        "4013034045057175611230"
        "4018087834845041231387857560"
        "4013034045057011231756"
        "2015056283484521231380"
        "2018087273474621231380"
        "1701134546716126236350"
        "2018087283484521231382785756"
        "2015056021241450"
        "2018087278575631241461670"
        "1701031241461670"
        "2018087031241461672785756"
        "2015056021231382834845"
        "2018087278575621231382834845"
        "1701046716126236351345"
        "2018087021231382734746"
        "3346360301001123"
        "4013034045057175601123"
        "4013034045057017561123"
        "4018087834845004123138785756"
        "2015056283484502123138"
        "2018087273474602123138"
        "2018087283484527857562123138"
        "1701134504671612623635"
        "2015056028348452123138"
        "2018087278575628348452123138"
        "2018087027347462123138"
        "1701013454671612623635"
        "2015056002124145"
        "2018087278575603124146167"
        "2018087027857563124146167"
        "1701003124146167"
        "3124146167170100"
        "2124145201505600"
        "3124146167201808727857560"
        "3124146167201808702785756"
        "4123138785756401808783484500"
        "1123334636030100"
        "1123401303404505717560"
        "1123401303404505701756"
        "4671612623635170113450"
        "2123138201505628348450"
        "2123138201808727347460"
        "2123138201808728348452785756"
        "4671612623635170101345"
        "2123138201505602834845"
        "2123138201808727857562834845"
        "2123138201808702734746"
        "2734746202303700"
        "1345356050252300"
        "2834845202303727857560"
        "2834845202303702785756"
        "1756370272474500"
        "0202404600"
        "0370272474517560"
        "0370272474501756"
        "2785756202303728348450"
        "0356050252313450"
        "0202303727347460"
        "0202303728348452785756"
        "2785756202303702834845"
        "0356050252301345"
        "0202303727857562834845"
        "0202303702734746"
        "2734746201808721231380"
        "2834845201505621231380"
        "1345170146716126236350"
        "2834845201808721231382785756"
        "1756401303404505711230"
        "0334636030111230"
        "0401808783484541231387857560"
        "0401303404505711231756"
        "2785756201808731241461670"
        "0201505621241450"
        "0170131241461670"
        "0201808731241461672785756"
        "2785756201808721231382834845"
        "0201505621231382834845"
        "0170146716126236351345"
        "0201808721231382734746"
        "2734746201808702123138"
        "2834845201505602123138"
        "2834845201808727857562123138"
        "1345170104671612623635"
        "1756401303404505701123"
        "0334636030101123"
        "0401303404505717561123"
        "0401808783484504123138785756"
        "2785756201808728348452123138"
        "0201505628348452123138"
        "0201808727347462123138"
        "0170113454671612623635"
        "2785756201808703124146167"
        "0201505602124145"
        "0201808727857563124146167"
        "0170103124146167"
        "3124146167017010"
        "3124146167278575620180870"
        "2124145020150560"
        "3124146167020180872785756"
        "4671612623635134517010"
        "2123138273474620180870"
        "2123138283484520150560"
        "2123138283484520180872785756"
        "4123138785756040180878348450"
        "1123175640130340450570"
        "1123033463603010"
        "1123040130340450571756"
        "4671612623635017011345"
        "2123138278575620180872834845"
        "2123138020150562834845"
        "2123138020180872734746"
        "2734746212313820180870"
        "1345467161262363517010"
        "2834845212313820150560"
        "2834845212313820180872785756"
        "2785756312414616720180870"
        "0312414616717010"
        "0212414520150560"
        "0312414616720180872785756"
        "1756112340130340450570"
        "0412313878575640180878348450"
        "0112333463603010"
        "0112340130340450571756"
        "2785756212313820180872834845"
        "0467161262363517011345"
        "0212313820150562834845"
        "0212313820180872734746"
        "2734746020230370"
        "2834845278575620230370"
        "1345035605025230"
        "2834845020230372785756"
        "2785756283484520230370"
        "0273474620230370"
        "0134535605025230"
        "0283484520230372785756"
        "1756037027247450"
        "0175637027247450"
        "0020240460"
        "0037027247451756"
        "2785756020230372834845"
        "0278575620230372834845"
        "0035605025231345"
        "0020230372734746"
        "2734746020180872123138"
        "2834845278575620180872123138"
        "2834845020150562123138"
        "1345017014671612623635"
        "2785756283484520180872123138"
        "0273474620180872123138"
        "0283484520150562123138"
        "0134517014671612623635"
        "1756040130340450571123"
        "0175640130340450571123"
        "0033463603011123"
        "0040180878348454123138785756"
        "2785756020180873124146167"
        "0278575620180873124146167"
        "0020150562124145"
        "0017013124146167"
        "3124146167001701"
        "3124146167278575602018087"
        "3124146167027857562018087"
        "2124145002015056"
        "4671612623635134501701"
        "2123138273474602018087"
        "2123138283484527857562018087"
        "2123138283484502015056"
        "4671612623635013451701"
        "2123138278575628348452018087"
        "2123138027347462018087"
        "2123138028348452015056"
        "4123138785756004018087834845"
        "1123175604013034045057"
        "1123017564013034045057"
        "1123003346360301"
        "2734746212313802018087"
        "1345467161262363501701"
        "2834845212313827857562018087"
        "2834845212313802015056"
        "2785756312414616702018087"
        "0312414616701701"
        "0312414616727857562018087"
        "0212414502015056"
        "2785756212313828348452018087"
        "0467161262363513451701"
        "0212313827347462018087"
        "0212313828348452015056"
        "1756112304013034045057"
        "0412313878575604018087834845"
        "0112317564013034045057"
        "0112303346360301"
        "2734746021231382018087"
        "2834845278575621231382018087"
        "1345046716126236351701"
        "2834845021231382015056"
        "2785756283484521231382018087"
        "0273474621231382018087"
        "0134546716126236351701"
        "0283484521231382015056"
        "2785756031241461672018087"
        "0278575631241461672018087"
        "0031241461671701"
        "0021241452015056"
        "1756011234013034045057"
        "0175611234013034045057"
        "0041231387857564018087834845"
        "0011233346360301"
        "2734746002023037"
        "2834845278575602023037"
        "2834845027857562023037"
        "1345003560502523"
        "2785756283484502023037"
        "0273474602023037"
        "0283484527857562023037"
        "0134503560502523"
        "2785756028348452023037"
        "0278575628348452023037"
        "0027347462023037"
        "0013453560502523"
        "1756003702724745"
        "0175603702724745"
        "0017563702724745"
        "0002024046";
    char const* quaternary_token = QUATERNARY_TABLE;

    par_msquares_quaternary_triangle_table = PAR_ALLOC(int*, 256);
    for (int i = 0; i < 256; i++) {
        int ntris = *quaternary_token++ - '0';
        for (int j = 0; j < ntris * 3; j++) {
            int pt = *quaternary_token++ - '0';
            assert(pt >= 0 && pt < 9);
        }
        ntris = *quaternary_token++ - '0';
        for (int j = 0; j < ntris * 3; j++) {
            int pt = *quaternary_token++ - '0';
            assert(pt >= 0 && pt < 9);
        }
        ntris = *quaternary_token++ - '0';
        for (int j = 0; j < ntris * 3; j++) {
            int pt = *quaternary_token++ - '0';
            assert(pt >= 0 && pt < 9);
        }
        ntris = *quaternary_token++ - '0';
        for (int j = 0; j < ntris * 3; j++) {
            int pt = *quaternary_token++ - '0';
            assert(pt >= 0 && pt < 9);
        }
    }
}

typedef struct {
    float const* data;
    float threshold;
    float lower_bound;
    float upper_bound;
    int width;
    int height;
} par_gray_context;

static int gray_inside(int location, void* contextptr)
{
    par_gray_context* context = (par_gray_context*) contextptr;
    return context->data[location] > context->threshold;
}

static int gray_multi_inside(int location, void* contextptr)
{
    par_gray_context* context = (par_gray_context*) contextptr;
    float val = context->data[location];
    float upper = context->upper_bound;
    float lower = context->lower_bound;
    return val >= lower && val < upper;
}

static float gray_height(float x, float y, void* contextptr)
{
    par_gray_context* context = (par_gray_context*) contextptr;
    int i = PAR_CLAMP(context->width * x, 0, context->width - 1);
    int j = PAR_CLAMP(context->height * y, 0, context->height - 1);
    return context->data[i + j * context->width];
}

typedef struct {
    par_byte const* data;
    par_byte color[4];
    int bpp;
    int width;
    int height;
} par_color_context;

static int color_inside(int location, void* contextptr)
{
    par_color_context* context = (par_color_context*) contextptr;
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
    par_color_context* context = (par_color_context*) contextptr;
    assert(context->bpp == 4);
    int i = PAR_CLAMP(context->width * x, 0, context->width - 1);
    int j = PAR_CLAMP(context->height * y, 0, context->height - 1);
    int k = i + j * context->width;
    return context->data[k * 4 + 3] / 255.0;
}

par_msquares_meshlist* par_msquares_color(par_byte const* data, int width,
    int height, int cellsize, uint32_t color, int bpp, int flags)
{
    par_color_context context;
    context.bpp = bpp;
    context.color[0] = (color >> 16) & 0xff;
    context.color[1] = (color >> 8) & 0xff;
    context.color[2] = (color & 0xff);
    context.color[3] = (color >> 24) & 0xff;
    context.data = data;
    context.width = width;
    context.height = height;
    return par_msquares_function(
        width, height, cellsize, flags, &context, color_inside, color_height);
}

par_msquares_meshlist* par_msquares_grayscale(float const* data, int width,
    int height, int cellsize, float threshold, int flags)
{
    par_gray_context context;
    context.width = width;
    context.height = height;
    context.data = data;
    context.threshold = threshold;
    return par_msquares_function(
        width, height, cellsize, flags, &context, gray_inside, gray_height);
}

par_msquares_meshlist* par_msquares_grayscale_multi(float const* data,
    int width, int height, int cellsize, float const* thresholds,
    int nthresholds, int flags)
{
    par_msquares_meshlist* mlists[2];
    mlists[0] = PAR_ALLOC(par_msquares_meshlist, 1);
    int connect = flags & PAR_MSQUARES_CONNECT;
    int snap = flags & PAR_MSQUARES_SNAP;
    int heights = flags & PAR_MSQUARES_HEIGHTS;
    if (!heights) {
        snap = connect = 0;
    }
    flags &= ~PAR_MSQUARES_INVERT;
    flags &= ~PAR_MSQUARES_DUAL;
    flags &= ~PAR_MSQUARES_CONNECT;
    flags &= ~PAR_MSQUARES_SNAP;
    par_gray_context context;
    context.width = width;
    context.height = height;
    context.data = data;
    context.lower_bound = -FLT_MAX;
    for (int i = 0; i <= nthresholds; i++) {
        int mergeconf = i > 0 ? connect : 0;
        if (i == nthresholds) {
            context.upper_bound = FLT_MAX;
            mergeconf |= snap;
        } else {
            context.upper_bound = thresholds[i];
        }
        mlists[1] = par_msquares_function(width, height, cellsize, flags,
            &context, gray_multi_inside, gray_height);
        mlists[0] = par_msquares_merge(mlists, 2, mergeconf);
        context.lower_bound = context.upper_bound;
        flags |= connect;
    }
    return mlists[0];
}

par_msquares_mesh const* par_msquares_get_mesh(
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

// Combine multiple meshlists by moving mesh pointers, and optionally apply
// a "snap" operation that assigns a single Z value across all verts in each
// mesh.  The Z value determined by the mesh's position in the final mesh list.
static par_msquares_meshlist* par_msquares_merge(par_msquares_meshlist** lists,
    int count, int snap)
{
    par_msquares_meshlist* merged = PAR_ALLOC(par_msquares_meshlist, 1);
    merged->nmeshes = 0;
    for (int i = 0; i < count; i++) {
        merged->nmeshes += lists[i]->nmeshes;
    }
    merged->meshes = PAR_ALLOC(par_msquares_mesh*, merged->nmeshes);
    par_msquares_mesh** pmesh = merged->meshes;
    for (int i = 0; i < count; i++) {
        par_msquares_meshlist* meshlist = lists[i];
        for (int j = 0; j < meshlist->nmeshes; j++) {
            *pmesh++ = meshlist->meshes[j];
        }
        free(meshlist);
    }
    if (!snap) {
        return merged;
    }
    pmesh = merged->meshes;
    float zmin = FLT_MAX;
    float zmax = -zmin;
    for (int i = 0; i < merged->nmeshes; i++, pmesh++) {
        float* pzed = (*pmesh)->points + 2;
        for (int j = 0; j < (*pmesh)->npoints; j++, pzed += 3) {
            zmin = PAR_MIN(*pzed, zmin);
            zmax = PAR_MAX(*pzed, zmax);
        }
    }
    float zextent = zmax - zmin;
    pmesh = merged->meshes;
    for (int i = 0; i < merged->nmeshes; i++, pmesh++) {
        float* pzed = (*pmesh)->points + 2;
        float zed = zmin + zextent * i / (merged->nmeshes - 1);
        for (int j = 0; j < (*pmesh)->npoints; j++, pzed += 3) {
            *pzed = zed;
        }
    }
    if (!(snap & PAR_MSQUARES_CONNECT)) {
        return merged;
    }
    for (int i = 1; i < merged->nmeshes; i++) {
        par_msquares_mesh* mesh = merged->meshes[i];

        // Find all extrusion points.  This is tightly coupled to the
        // tessellation code, which generates two "connector" triangles for each
        // extruded edge.  The first two verts of the second triangle are the
        // verts that need to be displaced.
        char* markers = PAR_ALLOC(char, mesh->npoints);
        int tri = mesh->ntriangles - mesh->nconntriangles;
        while (tri < mesh->ntriangles) {
            markers[mesh->triangles[tri * 3 + 3]] = 1;
            markers[mesh->triangles[tri * 3 + 4]] = 1;
            tri += 2;
        }

        // Displace all extrusion points down to the previous level.
        float zed = zmin + zextent * (i - 1) / (merged->nmeshes - 1);
        float* pzed = mesh->points + 2;
        for (int j = 0; j < mesh->npoints; j++, pzed += 3) {
            if (markers[j]) {
                *pzed = zed;
            }
        }
        free(markers);
    }
    return merged;
}

par_msquares_meshlist* par_msquares_function(int width, int height,
    int cellsize, int flags, void* context, par_msquares_inside_fn insidefn,
    par_msquares_height_fn heightfn)
{
    assert(width > 0 && width % cellsize == 0);
    assert(height > 0 && height % cellsize == 0);

    if (flags & PAR_MSQUARES_DUAL) {
        int connect = flags & PAR_MSQUARES_CONNECT;
        int snap = flags & PAR_MSQUARES_SNAP;
        int heights = flags & PAR_MSQUARES_HEIGHTS;
        if (!heights) {
            snap = connect = 0;
        }
        flags ^= PAR_MSQUARES_INVERT;
        flags &= ~PAR_MSQUARES_DUAL;
        flags &= ~PAR_MSQUARES_CONNECT;
        par_msquares_meshlist* m[2];
        m[0] = par_msquares_function(width, height, cellsize, flags,
            context, insidefn, heightfn);
        flags ^= PAR_MSQUARES_INVERT;
        if (connect) {
            flags |= PAR_MSQUARES_CONNECT;
        }
        m[1] = par_msquares_function(width, height, cellsize, flags,
            context, insidefn, heightfn);
        return par_msquares_merge(m, 2, snap | connect);
    }

    int invert = flags & PAR_MSQUARES_INVERT;

    // Create the two code tables if we haven't already.  These are tables of
    // fixed constants, so it's embarassing that we use dynamic memory
    // allocation for them.  However it's easy and it's one-time-only.

    if (!par_msquares_binary_point_table) {
        par_init_tables();
    }

    // Allocate the meshlist and the first mesh.

    par_msquares_meshlist* mlist = PAR_ALLOC(par_msquares_meshlist, 1);
    mlist->nmeshes = 1;
    mlist->meshes = PAR_ALLOC(par_msquares_mesh*, 1);
    mlist->meshes[0] = PAR_ALLOC(par_msquares_mesh, 1);
    par_msquares_mesh* mesh = mlist->meshes[0];
    mesh->dim = (flags & PAR_MSQUARES_HEIGHTS) ? 3 : 2;
    int ncols = width / cellsize;
    int nrows = height / cellsize;

    // Worst case is four triangles and six verts per cell, so allocate that
    // much.

    int maxtris = ncols * nrows * 4;
    int maxpts = ncols * nrows * 6;
    int maxedges = ncols * nrows * 2;

    // However, if we include extrusion triangles for boundary edges,
    // we need space for another 4 triangles and 4 points per cell.

    uint16_t* conntris = 0;
    int nconntris = 0;
    uint16_t* edgemap = 0;
    if (flags & PAR_MSQUARES_CONNECT) {
        conntris = PAR_ALLOC(uint16_t, maxedges * 6);
        maxtris +=  maxedges * 2;
        maxpts += maxedges * 2;
        edgemap = PAR_ALLOC(uint16_t, maxpts);
        for (int i = 0; i < maxpts; i++) {
            edgemap[i] = 0xffff;
        }
    }

    uint16_t* tris = PAR_ALLOC(uint16_t, maxtris * 3);
    int ntris = 0;
    float* pts = PAR_ALLOC(float, maxpts * mesh->dim);
    int npts = 0;

    // The "verts" x/y/z arrays are the 4 corners and 4 midpoints around the
    // square, in counter-clockwise order.  The origin of "triangle space" is at
    // the lower-left, although we expect the image data to be in raster order
    // (starts at top-left).

    float normalization = 1.0f / PAR_MAX(width, height);
    float normalized_cellsize = cellsize * normalization;
    int maxrow = (height - 1) * width;
    uint16_t* ptris = tris;
    uint16_t* pconntris = conntris;
    float* ppts = pts;
    float vertsx[8], vertsy[8];
    uint8_t* prevrowmasks = PAR_ALLOC(uint8_t, ncols);
    int* prevrowinds = PAR_ALLOC(int, ncols * 3);

    // If simplification is enabled, we need to track all 'F' cells and their
    // repsective triangle indices.
    uint8_t* simplification_codes = 0;
    uint16_t* simplification_tris = 0;
    uint8_t* simplification_ntris = 0;
    if (flags & PAR_MSQUARES_SIMPLIFY) {
        simplification_codes = PAR_ALLOC(uint8_t, nrows * ncols);
        simplification_tris = PAR_ALLOC(uint16_t, nrows * ncols);
        simplification_ntris = PAR_ALLOC(uint8_t, nrows * ncols);
    }

    // Do the march!
    for (int row = 0; row < nrows; row++) {
        vertsx[0] = vertsx[6] = vertsx[7] = 0;
        vertsx[1] = vertsx[5] = 0.5 * normalized_cellsize;
        vertsx[2] = vertsx[3] = vertsx[4] = normalized_cellsize;
        vertsy[0] = vertsy[1] = vertsy[2] = normalized_cellsize * (row + 1);
        vertsy[4] = vertsy[5] = vertsy[6] = normalized_cellsize * row;
        vertsy[3] = vertsy[7] = normalized_cellsize * (row + 0.5);

        int northi = row * cellsize * width;
        int southi = PAR_MIN(northi + cellsize * width, maxrow);
        int northwest = invert ^ insidefn(northi, context);
        int southwest = invert ^ insidefn(southi, context);
        int previnds[8] = {0};
        uint8_t prevmask = 0;

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

            int const* pointspec = par_msquares_binary_point_table[code];
            int ptspeclength = *pointspec++;
            int currinds[8] = {0};
            uint8_t mask = 0;
            uint8_t prevrowmask = prevrowmasks[col];
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

            int const* trianglespec = par_msquares_binary_triangle_table[code];
            int trispeclength = *trianglespec++;

            if (flags & PAR_MSQUARES_SIMPLIFY) {
                simplification_codes[ncols * row + col] = code;
                simplification_tris[ncols * row + col] = ntris;
                simplification_ntris[ncols * row + col] = trispeclength;
            }

            // Add triangles.
            while (trispeclength--) {
                int a = *trianglespec++;
                int b = *trianglespec++;
                int c = *trianglespec++;
                *ptris++ = currinds[c];
                *ptris++ = currinds[b];
                *ptris++ = currinds[a];
                ntris++;
            }

            // Create two extrusion triangles for each boundary edge.
            if (flags & PAR_MSQUARES_CONNECT) {
                trianglespec = par_msquares_binary_triangle_table[code];
                trispeclength = *trianglespec++;
                while (trispeclength--) {
                    int a = *trianglespec++;
                    int b = *trianglespec++;
                    int c = *trianglespec++;
                    int i = currinds[a];
                    int j = currinds[b];
                    int k = currinds[c];
                    int u = 0, v = 0, w = 0;
                    if ((a % 2) && (b % 2)) {
                        u = v = 1;
                    } else if ((a % 2) && (c % 2)) {
                        u = w = 1;
                    } else if ((b % 2) && (c % 2)) {
                        v = w = 1;
                    } else {
                        continue;
                    }
                    if (u && edgemap[i] == 0xffff) {
                        for (int d = 0; d < mesh->dim; d++) {
                            *ppts++ = pts[i * mesh->dim + d];
                        }
                        edgemap[i] = npts++;
                    }
                    if (v && edgemap[j] == 0xffff) {
                        for (int d = 0; d < mesh->dim; d++) {
                            *ppts++ = pts[j * mesh->dim + d];
                        }
                        edgemap[j] = npts++;
                    }
                    if (w && edgemap[k] == 0xffff) {
                        for (int d = 0; d < mesh->dim; d++) {
                            *ppts++ = pts[k * mesh->dim + d];
                        }
                        edgemap[k] = npts++;
                    }
                    if ((a % 2) && (b % 2)) {
                        *pconntris++ = i;
                        *pconntris++ = j;
                        *pconntris++ = edgemap[j];
                        *pconntris++ = edgemap[j];
                        *pconntris++ = edgemap[i];
                        *pconntris++ = i;
                    } else if ((a % 2) && (c % 2)) {
                        *pconntris++ = edgemap[k];
                        *pconntris++ = k;
                        *pconntris++ = i;
                        *pconntris++ = edgemap[i];
                        *pconntris++ = edgemap[k];
                        *pconntris++ = i;
                    } else if ((b % 2) && (c % 2)) {
                        *pconntris++ = j;
                        *pconntris++ = k;
                        *pconntris++ = edgemap[k];
                        *pconntris++ = edgemap[k];
                        *pconntris++ = edgemap[j];
                        *pconntris++ = j;
                    }
                    nconntris += 2;
                }
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
    free(edgemap);
    free(prevrowmasks);
    free(prevrowinds);

    // Perform quick-n-dirty simplification by iterating two rows at a time.
    // In no way does this create the simplest possible mesh, but at least it's
    // fast and easy.
    if (flags & PAR_MSQUARES_SIMPLIFY) {
        int in_run = 0, start_run;

        // First figure out how many triangles we can eliminate.
        int neliminated_triangles = 0;
        for (int row = 0; row < nrows - 1; row += 2) {
            for (int col = 0; col < ncols; col++) {
                int a = simplification_codes[ncols * row + col] == 0xf;
                int b = simplification_codes[ncols * row + col + ncols] == 0xf;
                if (a && b) {
                    if (!in_run) {
                        in_run = 1;
                        start_run = col;
                    }
                    continue;
                }
                if (in_run) {
                    in_run = 0;
                    int run_width = col - start_run;
                    neliminated_triangles += run_width * 4 - 2;
                }
            }
            if (in_run) {
                in_run = 0;
                int run_width = ncols - start_run;
                neliminated_triangles += run_width * 4 - 2;
            }
        }

        // Build a new index array cell-by-cell.  If any given cell is 'F' and
        // its neighbor to the south is also 'F', then it's part of a run.
        int nnewtris = ntris + nconntris - neliminated_triangles;
        uint16_t* newtris = PAR_ALLOC(uint16_t, nnewtris * 3);
        uint16_t* pnewtris = newtris;
        in_run = 0;
        for (int row = 0; row < nrows - 1; row += 2) {
            for (int col = 0; col < ncols; col++) {
                int cell = ncols * row + col;
                int south = cell + ncols;
                int a = simplification_codes[cell] == 0xf;
                int b = simplification_codes[south] == 0xf;
                if (a && b) {
                    if (!in_run) {
                        in_run = 1;
                        start_run = col;
                    }
                    continue;
                }
                if (in_run) {
                    in_run = 0;
                    int nw_cell = ncols * row + start_run;
                    int ne_cell = ncols * row + col - 1;
                    int sw_cell = nw_cell + ncols;
                    int se_cell = ne_cell + ncols;
                    int nw_tri = simplification_tris[nw_cell];
                    int ne_tri = simplification_tris[ne_cell];
                    int sw_tri = simplification_tris[sw_cell];
                    int se_tri = simplification_tris[se_cell];
                    int nw_corner = nw_tri * 3 + 4;
                    int ne_corner = ne_tri * 3 + 0;
                    int sw_corner = sw_tri * 3 + 2;
                    int se_corner = se_tri * 3 + 1;
                    *pnewtris++ = tris[se_corner];
                    *pnewtris++ = tris[sw_corner];
                    *pnewtris++ = tris[nw_corner];
                    *pnewtris++ = tris[nw_corner];
                    *pnewtris++ = tris[ne_corner];
                    *pnewtris++ = tris[se_corner];
                }
                int ncelltris = simplification_ntris[cell];
                int celltri = simplification_tris[cell];
                for (int t = 0; t < ncelltris; t++, celltri++) {
                    *pnewtris++ = tris[celltri * 3];
                    *pnewtris++ = tris[celltri * 3 + 1];
                    *pnewtris++ = tris[celltri * 3 + 2];
                }
                ncelltris = simplification_ntris[south];
                celltri = simplification_tris[south];
                for (int t = 0; t < ncelltris; t++, celltri++) {
                    *pnewtris++ = tris[celltri * 3];
                    *pnewtris++ = tris[celltri * 3 + 1];
                    *pnewtris++ = tris[celltri * 3 + 2];
                }
            }
            if (in_run) {
                in_run = 0;
                    int nw_cell = ncols * row + start_run;
                    int ne_cell = ncols * row + ncols - 1;
                    int sw_cell = nw_cell + ncols;
                    int se_cell = ne_cell + ncols;
                    int nw_tri = simplification_tris[nw_cell];
                    int ne_tri = simplification_tris[ne_cell];
                    int sw_tri = simplification_tris[sw_cell];
                    int se_tri = simplification_tris[se_cell];
                    int nw_corner = nw_tri * 3 + 4;
                    int ne_corner = ne_tri * 3 + 0;
                    int sw_corner = sw_tri * 3 + 2;
                    int se_corner = se_tri * 3 + 1;
                    *pnewtris++ = tris[se_corner];
                    *pnewtris++ = tris[sw_corner];
                    *pnewtris++ = tris[nw_corner];
                    *pnewtris++ = tris[nw_corner];
                    *pnewtris++ = tris[ne_corner];
                    *pnewtris++ = tris[se_corner];
            }
        }
        ptris = pnewtris;
        ntris -= neliminated_triangles;
        free(tris);
        tris = newtris;
        free(simplification_codes);
        free(simplification_tris);
        free(simplification_ntris);

        // Remove unreferenced points.
        char* markers = PAR_ALLOC(char, npts);
        ptris = tris;
        int newnpts = 0;
        for (int i = 0; i < ntris * 3; i++, ptris++) {
            if (!markers[*ptris]) {
                newnpts++;
                markers[*ptris] = 1;
            }
        }
        for (int i = 0; i < nconntris * 3; i++) {
            if (!markers[conntris[i]]) {
                newnpts++;
                markers[conntris[i]] = 1;
            }
        }
        float* newpts = PAR_ALLOC(float, newnpts * mesh->dim);
        uint16_t* mapping = PAR_ALLOC(uint16_t, (ntris + nconntris) * 3);
        ppts = pts;
        float* pnewpts = newpts;
        int j = 0;
        if (mesh->dim == 3) {
            for (int i = 0; i < npts; i++, ppts += 3) {
                if (markers[i]) {
                    *pnewpts++ = ppts[0];
                    *pnewpts++ = ppts[1];
                    *pnewpts++ = ppts[2];
                    mapping[i] = j++;
                }
            }
        } else {
            for (int i = 0; i < npts; i++, ppts += 2) {
                if (markers[i]) {
                    *pnewpts++ = ppts[0];
                    *pnewpts++ = ppts[1];
                    mapping[i] = j++;
                }
            }
        }
        free(pts);
        free(markers);
        pts = newpts;
        npts = newnpts;
        for (int i = 0; i < ntris * 3; i++) {
            tris[i] = mapping[tris[i]];
        }
        for (int i = 0; i < nconntris * 3; i++) {
            conntris[i] = mapping[conntris[i]];
        }
        free(mapping);
    }

    // Append all extrusion triangles to the main triangle array.
    // We need them to be last so that they form a contiguous sequence.
    pconntris = conntris;
    for (int i = 0; i < nconntris; i++) {
        *ptris++ = *pconntris++;
        *ptris++ = *pconntris++;
        *ptris++ = *pconntris++;
        ntris++;
    }
    free(conntris);

    // Final cleanup and return.
    assert(npts <= maxpts);
    assert(ntris <= maxtris);
    mesh->npoints = npts;
    mesh->points = pts;
    mesh->ntriangles = ntris;
    mesh->triangles = tris;
    mesh->nconntriangles = nconntris;
    return mlist;
}

#undef PAR_MIN
#undef PAR_MAX
#undef PAR_CLAMP
#undef PAR_ALLOC
#endif
