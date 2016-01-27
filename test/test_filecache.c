#include "describe.h"

#define PAR_FILECACHE_IMPLEMENTATION
#include "par_filecache.h"

#include <fcntl.h>
#include <unistd.h>

#if ENABLE_LZ4
#define PREFIX "build/tfc.lz4."
#else
#define PREFIX "build/tfc.unc."
#endif

int main()
{
    describe("the basics") {
        it("should be able to initialize and evict") {
            par_filecache_init(PREFIX, 1200);
            par_filecache_evict_all();
            par_filecache_evict_all();
        }
        it("should save and load small files") {
            char const* payload = "01234";
            par_filecache_save("small", (par_byte*) payload, 5, 0, 0);
            par_byte* received = 0;
            int nbytes = 0;
            int loaded = par_filecache_load("small", &received, &nbytes, 0, 0);
            assert_ok(loaded);
            assert_ok(received);
            assert_equal(nbytes, 5);
            assert_ok(!memcmp(received, payload, 5));
            free(received);
        }
        it("should save and load larger files") {
            char* payload = calloc(1024, 1);
            srand(1);
            for (int i = 0; i < 512; i++) {
                payload[i] = rand() % 256;
            }
            par_filecache_save("big", (par_byte*) payload, 1024, 0, 0);
            par_byte* received = 0;
            int nbytes = 0;
            int loaded = par_filecache_load("big", &received, &nbytes, 0, 0);
            assert_ok(loaded);
            assert_ok(received);
            assert_equal(nbytes, 1024);
            assert_ok(!memcmp(received, payload, 1024));
            free(received);
            free(payload);
        }
        it("should evict the oldest file when exceeding max size") {
            par_byte* received = 0;
            int nbytes = 0;
            char* payload = calloc(1024, 1);
            par_filecache_save("second", (par_byte*) payload, 1024, 0, 0);
            free(payload);
            int loaded = par_filecache_load("big", &received, &nbytes, 0, 0);
            #if ENABLE_LZ4
            assert_equal(loaded, 1);
            free(received);
            #else
            assert_equal(loaded, 0);
            #endif
        }
        it("returns false when not found") {
            par_byte* received = 0;
            int nbytes = 0;
            int loaded = par_filecache_load("bar", &received, &nbytes, 0, 0);
            assert_equal(loaded, 0);
        }
    }

    return assert_failures();
}
