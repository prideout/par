#include "par_asset.h"

#define PAR_MSQUARES_IMPLEMENTATION
#include "par_msquares.h"

int main(int argc, char* argv[])
{
    asset_init();

    int nbytes;
    par_byte* data;
    asset_get("msquares_color.png", &data, &nbytes);
    free(data);

    asset_get("msquares_island.1024.bin", &data, &nbytes);
    free(data);

    return 0;
}
