#include "par_asset.h"
#include "lodepng.h"

#define PAR_BLUENOISE_IMPLEMENTATION
#include "par_bluenoise.h"

static void test_bluenoise()
{
    int nbytes;
    par_byte* data;
    asset_get("msquares_color.png", &data, &nbytes);
    unsigned dims[2] = {0, 0};
    unsigned char* pixels;
    lodepng_decode_memory(&pixels, &dims[0], &dims[1], data, nbytes, LCT_RGBA,
        8);
    free(data);

    // DO STUFF

    free(pixels);
}

int main(int argc, char* argv[])
{
    asset_init();
    test_bluenoise();
    return 0;
}
