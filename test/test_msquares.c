#include "par_asset.h"
#include "lodepng.h"

#define PAR_MSQUARES_IMPLEMENTATION
#include "par_msquares.h"

#define CELLSIZE 32
#define IMGWIDTH 1024
#define IMGHEIGHT 1024
#define THRESHOLD 0.0f
#define OCEAN_COLOR 0x214562

static void test_color()
{
    int nbytes;
    par_byte* data;
    asset_get("msquares_color.png", &data, &nbytes);
    unsigned dims[2] = {0, 0};
    unsigned char* decoded;
    lodepng_decode_memory(&decoded, &dims[0], &dims[1], data, nbytes, LCT_RGBA,
        8);
    free(data);
    assert(dims[0] == IMGWIDTH);
    assert(dims[1] == IMGHEIGHT);
    int flags = PAR_MSQUARES_INVERT | PAR_MSQUARES_HEIGHTS;
    par_msquares_meshlist* mlist = par_msquares_from_color(decoded, IMGWIDTH,
        IMGHEIGHT, CELLSIZE, OCEAN_COLOR, 4, flags);
    free(decoded);
    // TODO generate obj here
    par_msquares_free(mlist);
}

static void test_grayscale()
{
    int nbytes;
    float* data;
    asset_get("msquares_island.1024.bin", (par_byte**) &data, &nbytes);
    int flags = PAR_MSQUARES_HEIGHTS;
    par_msquares_meshlist* mlist = par_msquares_from_grayscale(data,
        IMGWIDTH, IMGHEIGHT, CELLSIZE, THRESHOLD, flags);
    free(data);
    // TODO generate obj here
    par_msquares_free(mlist);
}

int main(int argc, char* argv[])
{
    asset_init();
    test_color();
    test_grayscale();
    return 0;
}
