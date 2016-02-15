#include "par_asset.h"
#include "lodepng.h"

#define PAR_BLUENOISE_IMPLEMENTATION
#include "par_bluenoise.h"

const int MAXPTS = 100000;
const int DENSITY = 200000;
const int RESOLUTION = 384;

#define CLAMP(x, min, max) ((x < min) ? min : ((x > max) ? max : x))

static void test_bluenoise()
{
    int nbytes;
    par_byte* data;

    // Download and decode the density image.
    asset_get("trillium.png", &data, &nbytes);
    unsigned dims[2] = {0, 0};
    unsigned char* pixels;
    lodepng_decode_memory(&pixels, &dims[0], &dims[1], data, nbytes, LCT_GREY,
        8);
    free(data);

    // Download the tileset and initialize the bluenoise context.
    nbytes = 0;
    asset_get("bluenoise.trimmed.bin", &data, &nbytes);
    par_bluenoise_context* ctx;
    ctx = par_bluenoise_from_buffer(data, nbytes, MAXPTS);
    free(data);

    // Copy the density image into the bluenoise context and free it.
    par_bluenoise_density_from_gray(ctx, pixels, dims[0], dims[1], 1);
    free(pixels);

    // Generate points.
    float left = -0.5;
    float bottom = -0.5;
    float right = 0.5;
    float top = 0.5;
    int npts;
    par_bluenoise_set_viewport(ctx, left, bottom, right, top);
    float* cpupts = par_bluenoise_generate(ctx, DENSITY, &npts);
    par_bluenoise_free(ctx);

    // Draw points.
    pixels = malloc(RESOLUTION * RESOLUTION);
    memset(pixels, 0xff, RESOLUTION * RESOLUTION);
    for (int i = 0; i < npts; i++) {
        float x = *cpupts++;
        float y = *cpupts++;
        cpupts++;
        int i = CLAMP(RESOLUTION * (x + 0.5), 0, RESOLUTION - 1);
        int j = CLAMP(RESOLUTION * (0.5 - y), 0, RESOLUTION - 1);
        pixels[i + j * RESOLUTION] = 0;
    }

    // Write out the image.
    lodepng_encode_file("build/bluenoise.png", pixels, RESOLUTION, RESOLUTION,
        LCT_GREY, 8);
    free(pixels);
}

int main(int argc, char* argv[])
{
    asset_init();
    test_bluenoise();
    return 0;
}
