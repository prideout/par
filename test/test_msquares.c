#include "par_asset.h"
#include "lodepng.h"

#define PAR_MSQUARES_IMPLEMENTATION
#include "par_msquares.h"

#define CELLSIZE 32
#define IMGWIDTH 1024
#define IMGHEIGHT 1024

int main(int argc, char* argv[])
{
    asset_init();

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
        IMGHEIGHT, CELLSIZE, 0x214562, 4, flags);
    free(decoded);

    // TODO generate obj here

    par_msquares_free(mlist);

    asset_get("msquares_island.1024.bin", &data, &nbytes);
    free(data);

    return 0;
}
