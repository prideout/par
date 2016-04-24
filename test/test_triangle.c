#include <stdio.h>

#define PAR_TRIANGLE_IMPLEMENTATION
#include "par_triangle.h"

static FILE* svg_begin(char const* path)
{
    FILE* svgfile = fopen(path, "wt");
    fputs(
        "<svg viewBox='-1.1 -1.1 10.2 10.2' width='500px' height='500px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n"
        "<g transform='translate(1.5 1.5) scale(0.3 0.3)'>", svgfile);
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
    int triangle, uint32_t color, float fill, float stroke)
{
    fprintf(svgfile, "<path\n"
        " fill='#%6.6x'\n"
        " stroke='#000000'\n"
        " stroke-width='0.05'\n"
        " fill-opacity='%f'\n"
        " stroke-opacity='%f'\n"
        " d='", color, fill, stroke);
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
    FILE* svgfile = svg_begin("build/triangle_simple.svg");
    par_triangle_mesh* mesh = par_triangle_mesh_create_cdt(path);
    for (int t = 0; t < mesh->ntriangles; t++) {
        svg_write_triangle(svgfile, mesh, t, 0xa0a0a0, 0.5, 0.5);
    }
    svg_write_path(svgfile, path, 0xff0000, 0.5, 0.5);
    par_triangle_path_free(path);
    par_triangle_mesh_free(mesh);
    svg_end(svgfile);
}

int main(int argc, char* argv[])
{
    test_simple();
    return 0;
}
