#include "describe.h"

#define PAR_BUBBLES_IMPLEMENTATION
#include "par_bubbles.h"

#define PAR_SHAPES_T uint32_t
#define PAR_SHAPES_IMPLEMENTATION
#include "par_shapes.h"

static par_bubbles_t* bubbles;

#define NRADIUSES 100
static double radiuses[NRADIUSES];

#define NNODES 252
static int hierarchy[NNODES] = {
      0,   0,   1,   2,   2,   2,   2,   1,   7,   7,   7,   7,
      7,   1,  13,   0,  15,  15,  15,  18,  18,  18,  18,  18,
     18,  18,  18,  18,  15,  15,  15,  15,  15,  15,  15,  15,
     15,   0,  37,  38,  38,  38,  38,  38,  37,  37,  37,  37,
     37,  37,   0,  50,  50,  50,  50,   0,  55,   0,  57,  57,
     57,  57,  57,  57,  57,  57,   0,  66,  66,  66,  66,  66,
     66,  66,  66,  66,  66,  66,  66,  66,  66,  66,  66,  66,
     66,  66,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,
     85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,
     85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  66,  66,
     66,  66,  66,  66,  66,  66,  66,  66,   0, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128,   0, 139, 139, 139, 139,
    139, 139, 139, 146, 146, 139, 139, 139, 139, 152, 152, 152,
    139, 139, 139, 158, 158, 158, 158, 139, 139, 139, 139, 139,
      0, 168, 169, 169, 169, 169, 169, 168, 175, 175, 175, 175,
    175, 175, 175, 175, 175, 175, 175, 168, 187, 187, 187, 187,
    187, 187, 193, 193, 193, 193, 187, 187, 187, 168, 201, 201,
    201, 201, 168, 206, 206, 206, 168, 210, 211, 211, 211, 210,
    215, 215, 215, 215, 215, 210, 221, 221, 221, 210, 210, 226,
    226, 226, 210, 230, 230, 230, 230, 230, 230, 230, 230, 230,
    230, 230, 230, 230, 230, 230, 210, 210, 210, 210, 210, 168,
};

static int get_depth(int* tree, int i)
{
    int d = 0;
    while (i) {
        i = tree[i];
        d++;
    }
    return d;
}

int main()
{
    for (int i = 0; i < NRADIUSES; i++) {
        radiuses[i] = 1 + rand() % 10;
    }

    describe("par_bubbles_pack") {

        it("should pass a simple smoke test") {
            bubbles = par_bubbles_pack(radiuses, NRADIUSES);
            assert_ok(bubbles);
            assert_equal(bubbles->count, (int) NRADIUSES);
            par_bubbles_export(bubbles, "build/test_bubbles_pack.svg");
            par_bubbles_free_result(bubbles);
        }

        it("should handle a small number of nodes") {
            bubbles = par_bubbles_pack(radiuses, 0);
            par_bubbles_export(bubbles, "build/test_bubbles_pack0.svg");
            par_bubbles_free_result(bubbles);
            bubbles = par_bubbles_pack(radiuses, 1);
            par_bubbles_export(bubbles, "build/test_bubbles_pack1.svg");
            par_bubbles_free_result(bubbles);
            bubbles = par_bubbles_pack(radiuses, 2);
            par_bubbles_export(bubbles, "build/test_bubbles_pack2.svg");
            par_bubbles_free_result(bubbles);
            bubbles = par_bubbles_pack(radiuses, 3);
            par_bubbles_export(bubbles, "build/test_bubbles_pack3.svg");
            par_bubbles_free_result(bubbles);
        }

        it("should work with small hierarchy") {
            static int hierarchy[10] = {
                0, 0, 0, 0, 1, 1, 2,
                5, 5, 5,
            };
            bubbles = par_bubbles_hpack_circle(hierarchy, 10, 100);
            par_bubbles_export(bubbles, "build/test_bubbles_hpack_circle1.svg");
            par_bubbles_free_result(bubbles);
        }

        it("can be exported to an SVG") {
            bubbles = par_bubbles_hpack_circle(hierarchy, NNODES, 100);
            par_bubbles_export(bubbles, "build/test_bubbles_hpack_circle2.svg");
            par_bubbles_free_result(bubbles);
        }

        it("can be exported to an OBJ") {

            const int nnodes = 2e3;
            const float zscale = 0.01;

            // First, generate a random tree.  Square the random parent pointers
            // to make the graph distribution a bit more interesting, and to
            // make it easier for humans to find deep portions of the tree.
            int* tree = malloc(sizeof(int) * nnodes);
            tree[0] = 0;
            for (int i = 1; i < nnodes; i++) {
                float a = (float) rand() / RAND_MAX;
                float b = (float) rand() / RAND_MAX;
                tree[i] = i * a * b;
            }

            // Perform circle packing.
            bubbles = par_bubbles_hpack_circle(tree, nnodes, 1.0);

            // Create template shape.
            float normal[3] = {0, 0, 1};
            float center[3] = {0, 0, 0};
            par_shapes_mesh* template = par_shapes_create_disk(1.0, 64,
                center, normal);

            // Merge each circle into the scene.
            par_shapes_mesh* scene = par_shapes_create_empty();
            double const* xyr = bubbles->xyr;
            par_shapes_mesh* clone = 0;
            for (int i = 0; i < bubbles->count; i++, xyr += 3) {
                float d = get_depth(tree, i);
                clone = par_shapes_clone(template, clone);
                par_shapes_scale(clone, xyr[2], xyr[2], 1.0);
                par_shapes_translate(clone, xyr[0], xyr[1], d * zscale);
                par_shapes_merge(scene, clone);
            }

            // Export the OBJ file.
            char const* const filename = "build/bubbles.obj";
            par_shapes_export(scene, filename);

            // Free memory.
            par_shapes_free_mesh(template);
            par_shapes_free_mesh(clone);
            par_shapes_free_mesh(scene);
            par_bubbles_free_result(bubbles);
            free(tree);
        }

    }

    describe("precision") {

        #define H1 10
        #define H2 20

        it("works well with relative coordinate systems") {
            bubbles = par_bubbles_hpack_local(hierarchy, NNODES);
            par_bubbles_export_local(bubbles, 0,
                "build/test_bubbles_hpack_local1.svg");
            par_bubbles_export_local(bubbles, 158,
                "build/test_bubbles_hpack_local2.svg");
            par_bubbles_free_result(bubbles);
        }

        it("is poor with deep nesting and non-local packing") {
            static int hierarchy[H2] = {0};
            for (int i = 1; i < H1; i++) {
                hierarchy[i] = i - 1;
            }
            for (int i = H1; i < H2; i++) {
                hierarchy[i] = H1 - 1;
            }
            bubbles = par_bubbles_hpack_circle(hierarchy, H2, 1);
            par_bubbles_export_local(bubbles, H1 - 1,
                "build/test_bubbles_hpack_local3.svg");
            par_bubbles_free_result(bubbles);
        }

        it("is great with deep nesting and local packing") {
            static int hierarchy[H2] = {0};
            for (int i = 1; i < H1; i++) {
                hierarchy[i] = i - 1;
            }
            for (int i = H1; i < H2; i++) {
                hierarchy[i] = H1 - 1;
            }
            bubbles = par_bubbles_hpack_local(hierarchy, H2);
            par_bubbles_export_local(bubbles, H1 - 1,
                "build/test_bubbles_hpack_local4.svg");
            par_bubbles_free_result(bubbles);
        }

        #undef H1
        #undef H2
    }

    describe("par_bubbles_find_local") {

        it("finds the smallest node that completely encloses a box") {
            bubbles = par_bubbles_hpack_local(hierarchy, NNODES);
            // This aabb encloses node 158, which means node 139 is the deepest
            // circle that encloses it.
            double cx = 40.804406 / 100.0f;
            double cy = -15.209295 / 100.0f;
            double r = 8.133015 / 100.0f;
            double aabb[4] = { cx - r, cy - r, cx + r, cy + r };
            int node = par_bubbles_find_local(bubbles, aabb, 0);
            par_bubbles_free_result(bubbles);
            assert_equal(node, 139);
        }

        it("ditto, but using a non-root coordinate system") {
            bubbles = par_bubbles_hpack_local(hierarchy, NNODES);
            par_bubbles_export_local(bubbles, 139,
                "build/test_bubbles_hpack_local5.svg");
            double cx = -0.511578, cy = 0.147831, r = 0.08;
            double aabb[4] = { cx - r, cy - r, cx + r, cy + r };
            int node = par_bubbles_find_local(bubbles, aabb, 139);
            par_bubbles_free_result(bubbles);
            assert_equal(node, 158);
        }

        it("is used by the pick_local function") {
            bubbles = par_bubbles_hpack_local(hierarchy, NNODES);
            double cx = -0.654258, cy = 0.065455;
            double minradius = 0;
            int node = par_bubbles_pick_local(bubbles, cx, cy, 139, minradius);
            assert_equal(node, 162);
            minradius = 0.095;
            node = par_bubbles_pick_local(bubbles, cx, cy, 139, minradius);
            assert_equal(node, 158);
            par_bubbles_free_result(bubbles);
        }
    }

    return assert_failures();
}
