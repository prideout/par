#include <stdio.h>

#define PAR_TRIANGLE_IMPLEMENTATION
#include "par_triangle.h"

static FILE* svg_begin(char const* path, float t, float s)
{
    FILE* svgfile = fopen(path, "wt");
    fprintf(svgfile,
        "<svg viewBox='-1.1 -1.1 10.2 10.2' width='500px' height='500px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g transform='translate(%f %f) scale(%f %f)'>", t, t, s, s);
    return svgfile;
}

static void svg_end(FILE* svgfile)
{
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
}

static void svg_write_path(FILE* svgfile, par_triangle_path const* polygon,
    uint32_t color, float fill, float stroke)
{
    fprintf(svgfile, "<path\n"
        " fill='#%6.6x'\n"
        " stroke='#%6.6x'\n"
        " stroke-width='0.005'\n"
        " fill-opacity='%f'\n"
        " stroke-opacity='%f'\n"
        " d='", color, color, fill, stroke);
    for (int c = 0; c < polygon->nloops; c++) {
        float const* loop = polygon->loops[c];
        int length = (int) polygon->lengths[c];
        uint16_t first = 1;
        for (uint16_t s = 0; s < length; s++) {
            fprintf(svgfile, "%c%f,%f", first ? 'M' : 'L', loop[s * 2],
                loop[s * 2 + 1]);
            first = 0;
        }
        fputs("Z", svgfile);
    }
    fputs("'\n/>", svgfile);
}

static void svg_write_triangle(FILE* svgfile, par_triangle_mesh const* mesh,
    int triangle, uint32_t color, float fill, float stroke, float thickness)
{
    fprintf(svgfile, "<path\n"
        " fill='#%6.6x'\n"
        " stroke='#000000'\n"
        " stroke-width='%f'\n"
        " fill-opacity='%f'\n"
        " stroke-opacity='%f'\n"
        " d='", color, thickness, fill, stroke);
    uint16_t const* tri = mesh->triangles + triangle * 3;
    float const* a = mesh->points + tri[0] * 2;
    float const* b = mesh->points + tri[1] * 2;
    float const* c = mesh->points + tri[2] * 2;
    fprintf(svgfile, "M%f,%fL%f,%fL%f,%fZ",
        a[0], a[1], b[0], b[1], c[0], c[1]);
    fputs("'\n/>", svgfile);
}

static void test_simple()
{
    const uint16_t lengths[2] = {5, 3};
    const float points[] = {
        0, 0, 6, 0, 6, 4, 3, 7, 0, 4, // CCW pentagon
        3, 3, 4, 1, 2, 1, // CW triangle
    };
    par_triangle_path* path = par_triangle_path_create(lengths, 2, points, 0);
    FILE* svgfile = svg_begin("build/triangle_simple.svg", 1.5, 0.3);
    par_triangle_mesh* mesh = par_triangle_mesh_create_cdt(path);
    for (int t = 0; t < mesh->ntriangles; t++) {
        svg_write_triangle(svgfile, mesh, t, 0xa0a0a0, 0.5, 0.5, 0.1);
    }
    svg_write_path(svgfile, path, 0xff0000, 0.5, 0.5);
    par_triangle_path_free(path);
    par_triangle_mesh_free(mesh);
    svg_end(svgfile);
}

static void test_concave()
{
    // Figure snarfed from:
    // https://github.com/prideout/polygon.js/blob/master/src/application.coffee

    const float points[] = {
        440, 502, 414, 438, 404, 355, 394, 298, 442, 278, 470, 265, 480, 252,
        489, 232, 508, 198, 467, 173, 437, 236, 395, 251, 361, 257, 324, 212,
        317, 170, 327, 150, 349, 125, 367, 82, 353, 56, 308, 22, 244, 40, 233,
        75, 258, 146, 278, 159, 299, 216, 282, 277, 228, 246, 168, 180, 159,
        167, 117, 207, 194, 249, 223, 277, 263, 304, 277, 385, 259, 406, 225,
        429, 217, 435, 159, 496, 293, 520, 284, 451, 315, 406, 323, 381, 351,
        391, 354, 421, 370, 458, 344, 487, 335, 535,
    };
    const uint16_t lengths[1] = {sizeof(points) / (2 * sizeof(float))};

    par_triangle_path* path = par_triangle_path_create(lengths, 1, points, 0);
    FILE* svgfile = svg_begin("build/triangle_concave.svg", 0.0, 0.015);
    par_triangle_mesh* mesh = par_triangle_mesh_create_cdt(path);
    for (int t = 0; t < mesh->ntriangles; t++) {
        svg_write_triangle(svgfile, mesh, t, 0xa0a0a0, 0.5, 0.5, 1.0);
    }
    svg_write_path(svgfile, path, 0xff0000, 0.5, 0.5);
    par_triangle_path_free(path);
    par_triangle_mesh_free(mesh);
    svg_end(svgfile);
}

int main(int argc, char* argv[])
{
    test_simple();
    test_concave();
    return 0;
}
