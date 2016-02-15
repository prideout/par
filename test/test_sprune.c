#define PAR_SPRUNE_IMPLEMENTATION
#include "par_sprune.h"
#include "describe.h"

static par_sprune_context* context;

static float boxes20[20 * 4] = {
    0.41, 0.45, 0.57, 0.57,
    0.7, 0.049, 0.74, 0.073,
    0.3, 0.55, 0.46, 0.61,
    0.31, 0.19, 0.62, 0.38,
    0.3, 0.1, 0.34, 0.15,
    0.23, 0.34, 0.23, 0.34,
    0.58, 0.22, 0.82, 0.61, // 6
    0.18, 0.71, 0.27, 0.72,
    0.69, 0.34, 0.81, 0.55,
    0.53, 0.1, 0.56, 0.17,
    0.33, 0.69, 0.39, 0.69,
    0.64, 0.46, 0.65, 0.46,
    0.35, 0.78, 0.4, 0.8,
    0.69, 0.63, 0.79, 0.68,
    0.7, 0.51, 0.77, 0.53,
    0.17, 0.79, 0.79, 0.81,
    0.66, 0.49, 0.67, 0.5, // 16
    0.36, 0.49, 0.41, 0.5,
    0.58, 0.53, 0.58, 0.54,
    0.25, 0.81, 0.51, 0.85,
};

static float boxes10[10 * 4] = {
    0.55, 0.38, 0.55, 0.38,
    0.5, 0.36, 0.51, 0.37,
    0.28, 0.2, 0.4, 0.24,
    0.39, 0.19, 0.4, 0.22,
    0.31, 0.29, 0.7, 0.43,
    0.019, 0.4, 0.02, 0.46,
    0.58, 0.44, 0.64, 0.75,
    0.39, 0.47, 0.48, 0.5,
    0.15, 0.91, 0.21, 0.93,
    0.62, 0.58, 0.63, 0.58,
};

static void export_svg(const char* filename, const char* color, float* boxes,
    int nboxes)
{
    FILE* svgfile = fopen(filename, "wt");
    fprintf(svgfile,
        "<svg viewBox='-.1 -.1 1.2 1.2' width='640px' height='640px' "
        "version='1.1' "
        "xmlns='http://www.w3.org/2000/svg'>\n");
    fprintf(svgfile,
        "<rect fill-opacity='0.1' stroke='none' fill='%s' x='0' y='0' "
        "width='100%%' height='100%%'/>\n", color);
    fprintf(svgfile,
        "<g stroke-width='0.001' stroke-opacity='0.5' stroke='black' "
        "fill-opacity='0.2' fill='#2A8BB6'>\n");
    for (int i = 0; i < nboxes; i++) {
        float minx = boxes[i * 4 + 0];
        float miny = boxes[i * 4 + 1];
        float maxx = boxes[i * 4 + 2];
        float maxy = boxes[i * 4 + 3];
        float w = maxx - minx;
        float h = maxy - miny;
        fprintf(svgfile,
            "<rect x='%f' y='%f' width='%f' height='%f'/>\n", minx, miny, w, h);
        fprintf(svgfile, "<text text-anchor='middle' stroke='none' fill-opacity='1.0' "
            "x='%f' y='%f' font-size='%f'>%d</text>\n",
            0.5 * (minx + maxx), 0.5 * (miny + maxy), 0.025, i);

    }
    fputs("</g>\n</svg>", svgfile);
    fclose(svgfile);
}

int main()
{
    describe("par_sprune_overlap") {

        it("should pass a simple smoke test") {
            context = par_sprune_overlap(boxes20, 20, 0);
            assert_ok(context);
            export_svg("build/test_sprune_overlap_20.svg", "#2AB68B",
                boxes20, 20);
            par_sprune_overlap(boxes10, 10, context);
            export_svg("build/test_sprune_overlap_10.svg", "#8B2AB6",
                boxes10, 10);
            par_sprune_free_context(context);
        }

    }

    return assert_failures();
}
