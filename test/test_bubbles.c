#include "describe.h"

#define PAR_BUBBLES_IMPLEMENTATION
#include "par_bubbles.h"

#define NRADIUSES 20

static double radiuses[NRADIUSES];
static par_bubbles_t* bubbles;

int main()
{
    for (int i = 0; i < NRADIUSES; i++) {
        radiuses[i] = 1 + rand() % 10;
    }

    describe("par_bubbles_pack") {

        it("should pass a simple smoke test") {
            bubbles = par_bubbles_pack(radiuses, NRADIUSES);
            assert_ok(bubbles);
            // assert_equal(bubbles->count, (int) NRADIUSES);
            par_bubbles_export(bubbles, "build/test_bubbles_pack.svg");
            par_bubbles_free_result(bubbles);
        }

        it("should handle a small number of nodes") {

            bubbles = par_bubbles_pack(radiuses, 0);
            assert_ok(bubbles);
            par_bubbles_export(bubbles, "build/test_bubbles_pack0.svg");
            par_bubbles_free_result(bubbles);

            bubbles = par_bubbles_pack(radiuses, 1);
            assert_ok(bubbles);
            par_bubbles_export(bubbles, "build/test_bubbles_pack1.svg");
            par_bubbles_free_result(bubbles);

            bubbles = par_bubbles_pack(radiuses, 2);
            assert_ok(bubbles);
            par_bubbles_export(bubbles, "build/test_bubbles_pack2.svg");
            par_bubbles_free_result(bubbles);

            bubbles = par_bubbles_pack(radiuses, 3);
            assert_ok(bubbles);
            par_bubbles_export(bubbles, "build/test_bubbles_pack3.svg");
            par_bubbles_free_result(bubbles);
        }
    }

    return assert_failures();
}
