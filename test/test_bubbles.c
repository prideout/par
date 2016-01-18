#include "describe.h"

#define PAR_BUBBLES_IMPLEMENTATION
#include "par_bubbles.h"

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

        it("should work with big hierarchy") {
            bubbles = par_bubbles_hpack_circle(hierarchy, NNODES, 100);
            par_bubbles_export(bubbles, "build/test_bubbles_hpack_circle2.svg");
            par_bubbles_free_result(bubbles);
        }

    }

    return assert_failures();
}
