#include "par_asset.h"
#include "lodepng.h"

#define PAR_MSQUARES_IMPLEMENTATION
#include "par_msquares.h"

#define CELLSIZE 32
#define IMGWIDTH 1024
#define IMGHEIGHT 1024
#define THRESHOLD 0.0f
#define OCEAN_COLOR 0x214562

static void test_multi()
{
    unsigned dims[2] = {0, 0};
    unsigned char* pixels;
    lodepng_decode_file(&pixels, &dims[0], &dims[1], "test/rgb.png", LCT_RGB,
        8);
    printf("%d x %d\n", dims[0], dims[1]);

    uint16_t* index;
    float* pt;
    par_msquares_meshlist* mlist;
    par_msquares_mesh const* mesh;
    FILE* objfile;

    mlist = par_msquares_color_multi(pixels, dims[0], dims[1], CELLSIZE / 2, 3,
        0);
    objfile = fopen("build/msquares_multi.obj", "wt");
    int offset = 1;
    for (int m = 0; m < par_msquares_get_count(mlist); m++) {
        mesh = par_msquares_get_mesh(mlist, m);
        pt = mesh->points;
        for (int i = 0; i < mesh->npoints; i++) {
            float z = mesh->dim > 2 ? pt[2] : 0;
            fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
            pt += mesh->dim;
        }
        index = mesh->triangles;
        for (int i = 0; i < mesh->ntriangles; i++) {
            int a = offset + *index++;
            int b = offset + *index++;
            int c = offset + *index++;
            fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
        offset += mesh->npoints;
    }
    fclose(objfile);
    par_msquares_free(mlist);

    free(pixels);
}

static void test_color()
{
    int nbytes;
    par_byte* data;
    asset_get("msquares_color.png", &data, &nbytes);
    unsigned dims[2] = {0, 0};
    unsigned char* pixels;
    lodepng_decode_memory(&pixels, &dims[0], &dims[1], data, nbytes, LCT_RGBA,
        8);
    free(data);
    assert(dims[0] == IMGWIDTH);
    assert(dims[1] == IMGHEIGHT);

    int flags, i;
    uint16_t* index;
    float* pt;
    par_msquares_meshlist* mlist;
    par_msquares_mesh const* mesh;
    FILE* objfile;

    // -----------------------------
    // msquares_color_default
    // -----------------------------
    flags = 0;
    mlist = par_msquares_color(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        OCEAN_COLOR, 4, flags);
    mesh = par_msquares_get_mesh(mlist, 0);
    objfile = fopen("build/msquares_color_default.obj", "wt");
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = 1 + *index++;
        int b = 1 + *index++;
        int c = 1 + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // -----------------------------
    // msquares_color_invert_heights
    // -----------------------------
    flags = PAR_MSQUARES_INVERT | PAR_MSQUARES_HEIGHTS;
    mlist = par_msquares_color(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        OCEAN_COLOR, 4, flags);
    mesh = par_msquares_get_mesh(mlist, 0);
    objfile = fopen("build/msquares_color_invert_heights.obj", "wt");
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = 1 + *index++;
        int b = 1 + *index++;
        int c = 1 + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // -----------------------------
    // msquares_color_dual_heights
    // -----------------------------
    flags = PAR_MSQUARES_DUAL | PAR_MSQUARES_HEIGHTS | PAR_MSQUARES_SNAP |
        PAR_MSQUARES_CONNECT | PAR_MSQUARES_SIMPLIFY;
    mlist = par_msquares_color(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        OCEAN_COLOR, 4, flags);
    mesh = par_msquares_get_mesh(mlist, 0);
    objfile = fopen("build/msquares_color_simplify.obj", "wt");
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    int offset = 1;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    offset = mesh->npoints + 1;
    mesh = par_msquares_get_mesh(mlist, 1);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    free(pixels);
}

static void test_grayscale()
{
    int nbytes;
    float* pixels;
    asset_get("msquares_island.1024.bin", (par_byte**) &pixels, &nbytes);

    int flags, i;
    uint16_t* index;
    float* pt;
    par_msquares_meshlist* mlist;
    par_msquares_mesh const* mesh;
    FILE* objfile;
    int offset = 1;

    // -----------------------------
    // msquares_gray_default
    // -----------------------------
    flags = 0;
    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        THRESHOLD, flags);
    assert(par_msquares_get_count(mlist) == 1);
    objfile = fopen("build/msquares_gray_default.obj", "wt");
    mesh = par_msquares_get_mesh(mlist, 0);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // -----------------------------
    // msquares_gray_invert
    // -----------------------------
    flags = PAR_MSQUARES_INVERT;
    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        THRESHOLD, flags);
    assert(par_msquares_get_count(mlist) == 1);
    objfile = fopen("build/msquares_gray_invert.obj", "wt");
    mesh = par_msquares_get_mesh(mlist, 0);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // -----------------------------
    // msquares_gray_dual
    // -----------------------------
    flags = PAR_MSQUARES_DUAL;
    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        THRESHOLD, flags);
    assert(par_msquares_get_count(mlist) == 2);
    objfile = fopen("build/msquares_gray_dual.obj", "wt");
    mesh = par_msquares_get_mesh(mlist, 0);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    offset = 1;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    offset = mesh->npoints + 1;
    mesh = par_msquares_get_mesh(mlist, 1);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // -----------------------------
    // msquares_gray_heights
    // -----------------------------
    flags = PAR_MSQUARES_HEIGHTS;
    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        THRESHOLD, flags);
    mesh = par_msquares_get_mesh(mlist, 0);
    objfile = fopen("build/msquares_gray_heights.obj", "wt");
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = 1 + *index++;
        int b = 1 + *index++;
        int c = 1 + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // ------------------------------------
    // msquares_gray_heights_dual_snap
    // ------------------------------------
    flags = PAR_MSQUARES_HEIGHTS | PAR_MSQUARES_DUAL | PAR_MSQUARES_SNAP;
    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        THRESHOLD, flags);
    objfile = fopen("build/msquares_gray_heights_dual_snap.obj", "wt");
    mesh = par_msquares_get_mesh(mlist, 0);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    offset = 1;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    offset = mesh->npoints + 1;
    mesh = par_msquares_get_mesh(mlist, 1);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // --------------------------------------------
    // msquares_gray_heights_dual_snap_connect
    // --------------------------------------------
    flags = PAR_MSQUARES_HEIGHTS | PAR_MSQUARES_DUAL | PAR_MSQUARES_SNAP |
         PAR_MSQUARES_CONNECT;
    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        THRESHOLD, flags);
    objfile = fopen("build/msquares_gray_heights_dual_snap_connect.obj", "wt");
    mesh = par_msquares_get_mesh(mlist, 0);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    offset = 1;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    offset = mesh->npoints + 1;
    mesh = par_msquares_get_mesh(mlist, 1);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // ------------------------------
    // msquares_gray_heights_simplify
    // ------------------------------
    flags = PAR_MSQUARES_SIMPLIFY;
    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        THRESHOLD, flags);
    objfile = fopen("build/msquares_gray_heights_simplify.obj", "wt");
    mesh = par_msquares_get_mesh(mlist, 0);
    pt = mesh->points;
    for (i = 0; i < mesh->npoints; i++) {
        float z = mesh->dim > 2 ? pt[2] : 0;
        fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
        pt += mesh->dim;
    }
    index = mesh->triangles;
    offset = 1;
    for (i = 0; i < mesh->ntriangles; i++) {
        int a = offset + *index++;
        int b = offset + *index++;
        int c = offset + *index++;
        fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // ------------------------------
    // msquares_gray_multi_simplify
    // ------------------------------
    float thresholds[] = {0.0, 0.1};
    mlist = par_msquares_grayscale_multi(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        thresholds, 2, PAR_MSQUARES_SIMPLIFY);
    objfile = fopen("build/msquares_gray_multi_simplify.obj", "wt");
    assert(par_msquares_get_count(mlist) == 3);
    for (int m = 0; m < 3; m++) {
        mesh = par_msquares_get_mesh(mlist, m);
        pt = mesh->points;
        for (i = 0; i < mesh->npoints; i++) {
            float z = mesh->dim > 2 ? pt[2] : 0;
            fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
            pt += mesh->dim;
        }
    }
    offset = 1;
    for (int m = 0; m < 3; m++) {
        mesh = par_msquares_get_mesh(mlist, m);
        index = mesh->triangles;
        for (i = 0; i < mesh->ntriangles; i++) {
            int a = offset + *index++;
            int b = offset + *index++;
            int c = offset + *index++;
            fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
        offset += mesh->npoints;
    }
    fclose(objfile);
    par_msquares_free(mlist);

    // ------------------------------
    // msquares_gray_multi_everything
    // ------------------------------
    mlist = par_msquares_grayscale_multi(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        thresholds, 2, PAR_MSQUARES_SIMPLIFY | PAR_MSQUARES_HEIGHTS |
        PAR_MSQUARES_SNAP | PAR_MSQUARES_CONNECT);
    objfile = fopen("build/msquares_gray_multi_everything.obj", "wt");
    assert(par_msquares_get_count(mlist) == 3);
    for (int m = 0; m < 3; m++) {
        mesh = par_msquares_get_mesh(mlist, m);
        pt = mesh->points;
        for (i = 0; i < mesh->npoints; i++) {
            float z = mesh->dim > 2 ? pt[2] : 0;
            fprintf(objfile, "v %f %f %f\n", pt[0], pt[1], z);
            pt += mesh->dim;
        }
    }
    offset = 1;
    for (int m = 0; m < 3; m++) {
        mesh = par_msquares_get_mesh(mlist, m);
        index = mesh->triangles;
        for (i = 0; i < mesh->ntriangles; i++) {
            int a = offset + *index++;
            int b = offset + *index++;
            int c = offset + *index++;
            fprintf(objfile, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
        offset += mesh->npoints;
    }
    fclose(objfile);
    par_msquares_free(mlist);

    free(pixels);
}

int main(int argc, char* argv[])
{
    asset_init();
    test_color();
    test_grayscale();
    test_multi();
    return 0;
}
