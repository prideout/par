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
    _baseurl = sdsnew("https://prideout.net/assets/");
    int length = wai_getExecutablePath(0, 0, 0);
    _exedir = sdsnewlen("", length);
    int dirlen;
    wai_getExecutablePath(_exedir, length, &dirlen);
    sdsrange(_exedir, 0, dirlen);
    sds prefix = sdscat(sdsdup(_exedir), "cache_");
    par_filecache_init(prefix, 10 * 1024 * 1024);
    sdsfree(prefix);
}

static void asset_get(char const* filename, par_byte** data, int* nbytes)
{
    int exists = par_filecache_load(filename, data, nbytes, 0, 0);
    if (exists) {
        return;
    }
    sds srcurl = sdscat(sdsdup(_baseurl), filename);
    printf("Downloading %s...\n", srcurl);
    int success = par_easycurl_to_memory(srcurl, data, nbytes);
    if (!success) {
        exit(1);
    }
    par_filecache_save(filename, *data, *nbytes, 0, 0);
    sdsfree(srcurl);
}
