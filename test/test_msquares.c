#include "par_asset.h"
#include "lodepng.h"

#define PAR_MSQUARES_IMPLEMENTATION
#include "par_msquares.h"

#define CELLSIZE 32
#define IMGWIDTH 1024
#define IMGHEIGHT 1024
#define THRESHOLD 0.0f
#define OCEAN_COLOR 0x214562

static void svg_begin(FILE* svgfile)
{
    fputs(
        "<svg viewBox='-.1 -.1 1.2 1.2' width='500px' height='500px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g transform='translate(0 1) scale(1 -1)'>", svgfile);
}

static void svg_write_path(FILE* svgfile, par_msquares_mesh const* mesh,
    uint32_t color, float fill, float stroke)
{
    par_msquares_boundary* polygon = par_msquares_extract_boundary(mesh);
    fprintf(svgfile, "<path\n"
        " fill='#%6.6x'\n"
        " stroke='#%6.6x'\n"
        " stroke-width='0.005'\n"
        " fill-opacity='%f'\n"
        " stroke-opacity='%f'\n"
        " d='", color, color, fill, stroke);
    for (int c = 0; c < polygon->nchains; c++) {
        float const* chain = polygon->chains[c];
        int length = (int) polygon->lengths[c];
        uint16_t first = 1;
        for (uint16_t s = 0; s < length; s++) {
            fprintf(svgfile, "%c%f,%f", first ? 'M' : 'L', chain[s * 2],
                chain[s * 2 + 1]);
            first = 0;
        }
        fputs("Z", svgfile);
    }
    fputs("'\n/>", svgfile);
    par_msquares_free_boundary(polygon);
}

static void test_multi()
{
    unsigned dims[2] = {0, 0};
    int offset;
    unsigned char* pixels;
    uint16_t* index;
    float* pt;
    par_msquares_meshlist* mlist;
    par_msquares_mesh const* mesh;
    FILE* objfile, *svgfile;

    lodepng_decode_file(&pixels, &dims[0], &dims[1], "test/tverts.png", LCT_RGB,
        8);
    mlist = par_msquares_color_multi(pixels, dims[0], dims[1], 12, 3,
        PAR_MSQUARES_CLEAN);
    objfile = fopen("build/msquares_multi_tverts.obj", "wt");
    svgfile = fopen("build/msquares_multi_tverts.svg", "wt");
    offset = 1;
    svg_begin(svgfile);
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
        svg_write_path(svgfile, mesh, mesh->color & 0xffffff, 0.5, 0.5);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
    fclose(objfile);
    par_msquares_free(mlist);
    free(pixels);

    lodepng_decode_file(&pixels, &dims[0], &dims[1], "test/rgb.png", LCT_RGB,
        8);
    mlist = par_msquares_color_multi(pixels, dims[0], dims[1], CELLSIZE / 2, 3,
        PAR_MSQUARES_CLEAN);
    objfile = fopen("build/msquares_multi_rgb.obj", "wt");
    svgfile = fopen("build/msquares_multi_rgb.svg", "wt");
    offset = 1;
    svg_begin(svgfile);
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
        svg_write_path(svgfile, mesh, mesh->color, 0.5, 0.5);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
    fclose(objfile);
    par_msquares_free(mlist);
    free(pixels);

    lodepng_decode_file(&pixels, &dims[0], &dims[1], "test/rgba.png", LCT_RGBA,
        8);
    mlist = par_msquares_color_multi(pixels, dims[0], dims[1], CELLSIZE / 2, 4,
        PAR_MSQUARES_HEIGHTS | PAR_MSQUARES_CONNECT | PAR_MSQUARES_CLEAN);
    objfile = fopen("build/msquares_multi_rgba.obj", "wt");
    svgfile = fopen("build/msquares_multi_rgba.svg", "wt");
    offset = 1;
    svg_begin(svgfile);
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
        svg_write_path(svgfile, mesh, mesh->color & 0xffffff, 0.5, 0.5);
    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
    fclose(objfile);
    par_msquares_free(mlist);
    free(pixels);

    int nbytes;
    par_byte* data;
    asset_get("msquares_color.png", &data, &nbytes);
    lodepng_decode_memory(&pixels, &dims[0], &dims[1], data, nbytes, LCT_RGBA,
        8);
    free(data);
    mlist = par_msquares_color_multi(pixels, dims[0], dims[1], CELLSIZE, 4,
        PAR_MSQUARES_HEIGHTS | PAR_MSQUARES_CONNECT | PAR_MSQUARES_SIMPLIFY);
    objfile = fopen("build/msquares_multi_diagram.obj", "wt");
    offset = 1;
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
    FILE* objfile, *svgfile;
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
    svgfile = fopen("build/msquares_gray_dual.svg", "wt");
    svg_begin(svgfile);
    mesh = par_msquares_get_mesh(mlist, 0);
    svg_write_path(svgfile, mesh, 0x050b0, 0.5, 1.0);
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
    svg_write_path(svgfile, mesh, 0x10e050, 1.0, 0.0);
    svg_write_path(svgfile, mesh, 0, 0.0, 1.0);
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

    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        -0.02, flags);
    mesh = par_msquares_get_mesh(mlist, 0);
    svg_write_path(svgfile, mesh, 0xffffff, 0.0, 1.0);
    par_msquares_free(mlist);

    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        -0.06, flags);
    mesh = par_msquares_get_mesh(mlist, 0);
    svg_write_path(svgfile, mesh, 0xffffff, 0.0, 0.5);
    par_msquares_free(mlist);

    mlist = par_msquares_grayscale(pixels, IMGWIDTH, IMGHEIGHT, CELLSIZE,
        -0.1, flags);
    mesh = par_msquares_get_mesh(mlist, 0);
    svg_write_path(svgfile, mesh, 0xffffff, 0.0, 0.25);
    par_msquares_free(mlist);

    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);

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
