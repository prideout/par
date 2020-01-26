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

FILE* global_file;

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
                    "void main() { ... }\n\n");
            assert_str_equal(parsb_get_blocks(blocks, "common"),
                    "uniform vec4 resolution;\nuniform vec4 color;\n");
            assert_str_equal(parsb_get_blocks(blocks, "spooky"),
                    "-- hello\n!! world\n ");
        }

        it("can glue blocks") {
            assert_str_equal(parsb_get_blocks(blocks, "spooky common"),
                    "-- hello\n!! world\n uniform vec4 resolution;\nuniform vec4 color;\n");
        }

        it("can add a named block") {
            parsb_add_block(blocks, "prefix", "13");
            assert_str_equal(parsb_get_blocks(blocks, "prefix spooky common"),
                    "13-- hello\n!! world\n uniform vec4 resolution;\nuniform vec4 color;\n");
        }

        it("can replace a named block") {
            parsb_add_block(blocks, "prefix", "11");
            assert_str_equal(parsb_get_blocks(blocks, "prefix spooky common"),
                    "11-- hello\n!! world\n uniform vec4 resolution;\nuniform vec4 color;\n");
        }

        it("add and replace multiple blocks") {
            const char newblocks[] = R"(
--- prefix
12
--- great
goodbye)";
            parsb_add_blocks(blocks, newblocks, strlen(newblocks));
            assert_str_equal(parsb_get_blocks(blocks, "prefix spooky"),
                    "12\n-- hello\n!! world\n ");
            assert_str_equal(parsb_get_blocks(blocks, "great prefix"), "goodbye12\n");
        }

        it("can export the database") {
            global_file = fopen("test0.glsl", "w");
            parsb_write_blocks(blocks, [](const char* line, void* user) {
                fprintf(global_file, "%s\n", line);
            }, nullptr);
            fclose(global_file);
            parsb_destroy_context(blocks);
        }

        it("can import the database") {
            blocks = parsb_create_context((parsb_options){});
            assert_equal((int) sizeof(test_string), (int) strlen(test_string) + 1);
            parsb_add_blocks_from_file(blocks, "test0.glsl");
        }

        it("can re-export the database") {
            global_file = fopen("test1.glsl", "w");
            parsb_write_blocks(blocks, [](const char* line, void* user) {
                fprintf(global_file, "%s\n", line);
            }, nullptr);
            fclose(global_file);
            parsb_destroy_context(blocks);
        }
    }

    return assert_failures();
}
