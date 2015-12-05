#include "whereami.h"
#include "sds.h"

#define PAR_EASYCURL_IMPLEMENTATION
#include "par_easycurl.h"

#define PAR_FILECACHE_IMPLEMENTATION
#include "par_filecache.h"

static sds _exedir = 0;
static sds _baseurl = 0;

static void asset_init()
{
    par_easycurl_init(0);
    _baseurl = sdsnew("http://github.prideout.net/assets/");
    int length = wai_getExecutablePath(0, 0, 0);
    _exedir = sdsnewlen("", length);
    int dirlen;
    wai_getExecutablePath(_exedir, length, &dirlen);
    sdsrange(_exedir, 0, dirlen);
    par_filecache_init(_exedir, 10 * 1024 * 1024);
}

static void _fetch(char const* filename, char const* targetfolder)
{
    sds srcurl = sdscat(sdsdup(_baseurl), filename);
    sds dstpath = sdscat(sdsnew(targetfolder), filename);
    printf("Downloading %s...\n", srcurl);
    int success = par_easycurl_to_file(srcurl, dstpath);
    if (success) {
        // par_filecache_save(filename, payload, payloadsize, 0, 0);
    } else {
        exit(1);
    }
    sdsfree(srcurl);
    sdsfree(dstpath);
}

int main(int argc, char* argv[])
{
    asset_init();
    _fetch("msquares_color.png", _exedir);
    _fetch("msquares_island.1024.bin", _exedir);

    return 0;
}
