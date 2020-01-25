#define PAR_STRING_BLOCKS_IMPLEMENTATION
#include "par_string_blocks.h"

extern "C" {
#include "describe.h"
}

const char test_string[] = R"(
--- my_shader
void main() { ... }

--- common
uniform vec4 resolution;
uniform vec4 color;




--- spooky
-- hello
!! world


)";

int main()
{
    parsb_context* blocks;

    describe("simple") {

        it("can create a context") {
            blocks = parsb_create_context((parsb_options){});
            assert_equal((int) sizeof(test_string), (int) strlen(test_string) + 1);
            parsb_add_blocks(blocks, test_string, strlen(test_string));
        }

        it("can extract single blocks") {
            assert_str_equal(parsb_get_blocks(blocks, "my_shader"),
                    "void main() { ... }");
            assert_str_equal(parsb_get_blocks(blocks, "common"),
                    "uniform vec4 resolution;\nuniform vec4 color;");
            assert_str_equal(parsb_get_blocks(blocks, "spooky"),
                    "-- hello\n!! world");
        }

        it("can glue blocks") {
            assert_str_equal(parsb_get_blocks(blocks, "spooky common"),
                    "-- hello\n!! worlduniform vec4 resolution;\nuniform vec4 color;");
        }

        it("can add a named block") {
            parsb_add_block(blocks, "prefix", "13");
            assert_str_equal(parsb_get_blocks(blocks, "prefix spooky common"),
                    "13-- hello\n!! worlduniform vec4 resolution;\nuniform vec4 color;");
        }

        it("can replace a named block") {
            parsb_add_block(blocks, "prefix", "11");
            assert_str_equal(parsb_get_blocks(blocks, "prefix spooky common"),
                    "11-- hello\n!! worlduniform vec4 resolution;\nuniform vec4 color;");
        }

        it("add and replace multiple blocks") {
            const char newblocks[] = R"(
--- prefix
12
--- great
goodbye
            )";
            parsb_add_blocks(blocks, newblocks, strlen(newblocks));
            assert_str_equal(parsb_get_blocks(blocks, "prefix spooky"),
                    "12\n-- hello\n!! world");
            assert_str_equal(parsb_get_blocks(blocks, "great prefix"), "goodbye12\n");
        }

        it("can destroy the context") {
            parsb_destroy_context(blocks);
        }
    }

    return assert_failures();
}
