// FILECACHE :: https://github.com/prideout/par
// Simple file-based LRU cache for blobs with content-addressable names.
//
// Each cached item is stored on disk as "{PREFIX}{NAME}", where {PREFIX}
// is passed in when initializing the cache.  You'll probably want to specify
// a folder path for your prefix, including the trailing slash.
//
// Each item is divided into a payload (arbitrary size) and a header (a priori
// size). The structure of the payload and header are completely up to you. The
// list of items is stored in a text file at "{PREFIX}table", which contains a
// list of names, timestamps, and byte counts.  This table is loaded only once,
// but is saved every time the client fetches a blob from the cache, so that
// the most-recently-accessed timestamps are always up to date, even when
// your application doesn't close gracefully.
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

// Set this to zero if you wish to avoid LZ4 compression.  I recommend using
// it though, because it's very fast and it's a two-file library.
#ifndef ENABLE_LZ4
#define ENABLE_LZ4 0
#endif

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

typedef unsigned char par_byte;

// Initialize the filecache using the given prefix (usually a folder path with
// a trailing slash) and the given maximum byte count.  If items already exist
// in the cache when this is called, they are not evicted.  Cached items are
// meant to persist from run to run.
void par_filecache_init(const char* prefix, int maxsize);

// Save a blob to the cache using the given unique name.  If adding the blob
// would cause the cache to exceed maxsize, the least-recently-used item is
// evicted at this time.
void par_filecache_save(const char* name, par_byte* payload, int payloadsize,
    par_byte* header, int headersize);

// Check if the given blob is in the cache; if not, return 0.  If so, return 1
// and allocate new memory for payload.  The caller should free the payload.
// The header is preallocated so the caller needs to know its size beforehand.
int par_filecache_load(const char* name, par_byte** payload, int* payloadsize,
    par_byte* header, int headersize);

// Remove all items from the cache.
void par_filecache_evict_all();

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#include <limits.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <libgen.h>
#include <time.h>
#include <sys/stat.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#if ENABLE_LZ4
#include "lz4.h"
#endif

static char * _par_strdup(const char *s)
{
    if (s) {
        size_t l = strlen(s);
        char *s2 = malloc(l + 1);
        if (s2) {
            strcpy(s2, s);
        }
        return s2;
    }
    return 0;
}

#define MAX_ENTRIES 64
#define MIN(a, b) (a < b ? a : b)

typedef struct {
    time_t last_used_timestamp;
    uint64_t hashed_name;
    const char* name;
    int nbytes;
} filecache_entry_t;

typedef struct {
    filecache_entry_t entries[MAX_ENTRIES];
    int nentries;
    int totalbytes;
} filecache_table_t;

static void _update_table(const char* item_name, int item_size);
static void _append_table(const char* item_name, int item_size);
static void _read_or_create_tablefile();
static void _save_tablefile();
static void _evict_lru();
static uint64_t _hash(const char* name);

static char _fileprefix[PATH_MAX] = "./_cache.";
static char _tablepath[PATH_MAX] = "./_cache.table";
static int _maxtotalbytes = 1024 * 1024 * 16;
static filecache_table_t* _table = 0;

void par_filecache_init(const char* prefix, int maxsize)
{
    size_t len = strlen(prefix);
    assert(len + 1 < PATH_MAX && "Cache prefix is too long");
    strncpy(_fileprefix, prefix, len + 1);
    strcpy(_tablepath, _fileprefix);
    strcat(_tablepath, "table");
    _maxtotalbytes = maxsize;
}

#if IOS_EXAMPLE
NSString* getPrefix()
{
    NSString* cachesFolder = [NSSearchPathForDirectoriesInDomains(
        NSCachesDirectory, NSUserDomainMask, YES) firstObject];
    NSError* error = nil;
    if (![[NSFileManager defaultManager] createDirectoryAtPath : cachesFolder
        withIntermediateDirectories : YES
        attributes : nil
        error : &error]) {
        NSLog(@ "MGMPlatformGetCachesFolder error: %@", error);
        return nil;
    }
    return [cachesFolder stringByAppendingString : @ "/_cache."];
}

#endif

int par_filecache_load(const char* name, par_byte** payload, int* payloadsize,
    par_byte* header, int headersize)
{
    char qualified[PATH_MAX];
    size_t len = strlen(name);
    if (len == 0) {
        return 0;
    }
    assert(len + strlen(_fileprefix) < PATH_MAX);
    strcpy(qualified, _fileprefix);
    strcat(qualified, name);
    if (access(qualified, F_OK) == -1) {
        return 0;
    }
    FILE* cachefile = fopen(qualified, "rb");
    assert(cachefile && "Unable to open cache file for reading");
    fseek(cachefile, 0, SEEK_END);
    long fsize = ftell(cachefile);
    fseek(cachefile, 0, SEEK_SET);
    if (headersize > 0) {
        fread(header, headersize, 1, cachefile);
    }
    int32_t dnbytes;
    long cnbytes = fsize - headersize - sizeof(dnbytes);
    fread(&dnbytes, 1, sizeof(dnbytes), cachefile);
    char* cbuff = malloc(cnbytes);
    fread(cbuff, 1, cnbytes, cachefile);
#if ENABLE_LZ4
    char* dbuff = malloc(dnbytes);
    LZ4_decompress_safe(cbuff, dbuff, (int) cnbytes, dnbytes);
    free(cbuff);
#else
    char* dbuff = cbuff;
    dnbytes = (int32_t) cnbytes;
#endif
    fclose(cachefile);
    *payload = (par_byte*) dbuff;
    *payloadsize = dnbytes;
    _update_table(name, (int) cnbytes);
    return 1;
}

void par_filecache_save(const char* name, par_byte* payload, int payloadsize,
    par_byte* header, int headersize)
{
    char qualified[PATH_MAX];
    size_t len = strlen(name);
    if (len == 0) {
        return;
    }
    assert(len + strlen(_fileprefix) < PATH_MAX);
    strcpy(qualified, _fileprefix);
    strcat(qualified, name);
    FILE* cachefile = fopen(qualified, "wb");
    assert(cachefile && "Unable to open cache file for writing");
    if (headersize > 0) {
        fwrite(header, 1, headersize, cachefile);
    }
    int csize = 0;
    if (payloadsize > 0) {
        int32_t nbytes = payloadsize;
        fwrite(&nbytes, 1, sizeof(nbytes), cachefile);
#if ENABLE_LZ4
        int maxsize = LZ4_compressBound(nbytes);
        char* dst = malloc(maxsize);
        const char* src = (const char*) payload;
        assert(nbytes < LZ4_MAX_INPUT_SIZE);
        csize = LZ4_compress_default(src, dst, nbytes, maxsize);
        fwrite(dst, 1, csize, cachefile);
        free(dst);
#else
        csize = payloadsize;
        size_t actual = fwrite(payload, 1, csize, cachefile);
        if (actual < csize) {
            fclose(cachefile);
            remove(qualified);
            printf("Unable to save %s to cache (%d bytes)\n", name, csize);
            return;
        }
#endif
    }
    fclose(cachefile);
    _update_table(name, csize + headersize);
}

void par_filecache_evict_all()
{
    printf("Evicting all.\n");
    char qualified[PATH_MAX];
    if (!_table) {
        _read_or_create_tablefile();
    }
    filecache_entry_t* entry = _table->entries;
    for (int i = 0; i < _table->nentries; i++, entry++) {
        strcpy(qualified, _fileprefix);
        strcat(qualified, entry->name);
        printf("Evicting %s\n", qualified);
        remove(qualified);
    }
    _table->nentries = 0;
    _table->totalbytes = 0;
    remove(_tablepath);
}

// Adds the given item to the table and evicts the LRU items if the total cache
// size exceeds the specified maxsize.
static void _append_table(const char* item_name, int item_size)
{
    time_t now = time(0);
    if (!_table) {
        _read_or_create_tablefile();
    }
    int64_t hashed_name = _hash(item_name);
    int total = _table->totalbytes + item_size;
    while (_table->nentries >= MAX_ENTRIES || total > _maxtotalbytes) {
        assert(_table->nentries > 0 && "Cache size is too small.");
        _evict_lru();
        total = _table->totalbytes + item_size;
    }
    _table->totalbytes = total;
    filecache_entry_t* entry = &_table->entries[_table->nentries++];
    entry->last_used_timestamp = now;
    entry->hashed_name = hashed_name;
    entry->name = _par_strdup(item_name);
    entry->nbytes = item_size;
    _save_tablefile();
}

// Updates the timestamp associated with the given item.
static void _update_table(const char* item_name, int item_size)
{
    time_t now = time(0);
    if (!_table) {
        _read_or_create_tablefile();
    }
    int64_t hashed_name = _hash(item_name);
    filecache_entry_t* entry = _table->entries;
    int i;
    for (i = 0; i < _table->nentries; i++, entry++) {
        if (entry->hashed_name == hashed_name) {
            break;
        }
    }
    if (i >= _table->nentries) {
        _append_table(item_name, item_size);
        return;
    }
    entry->last_used_timestamp = now;
    _save_tablefile();
}

static void _read_or_create_tablefile()
{
    _table = calloc(sizeof(filecache_table_t), 1);
    FILE* fhandle = fopen(_tablepath, "r");
    if (!fhandle) {
        fhandle = fopen(_tablepath, "w");
        if (!fhandle) {
            mkdir(dirname(_tablepath), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            fhandle = fopen(_tablepath, "w");
        }
        assert(fhandle && "Unable to create filecache info file.");
    } else {
        filecache_entry_t entry;
        char name[PATH_MAX];
        while (1) {
            int nargs = fscanf(fhandle, "%ld %d %s", &entry.last_used_timestamp,
                    &entry.nbytes, name);
            if (nargs != 3) {
                break;
            }
            entry.name = _par_strdup(name);
            entry.hashed_name = _hash(entry.name);
            _table->entries[_table->nentries++] = entry;
            _table->totalbytes += entry.nbytes;
        }
    }
    fclose(fhandle);
}

static void _save_tablefile()
{
    FILE* fhandle = fopen(_tablepath, "w");
    assert(fhandle && "Unable to create filecache info file.");
    filecache_entry_t* entry = _table->entries;
    for (int i = 0; i < _table->nentries; i++, entry++) {
        fprintf(fhandle, "%ld %d %s\n", entry->last_used_timestamp,
            entry->nbytes, entry->name);
    }
    fclose(fhandle);
}

static void _evict_lru()
{
    const int64_t never_evict = _hash("version");
    int oldest_index = -1;
    time_t oldest_time = LONG_MAX;
    filecache_entry_t* entry = _table->entries;
    for (int i = 0; i < _table->nentries; i++, entry++) {
        if (entry->hashed_name == never_evict) {
            continue;
        }
        if (entry->last_used_timestamp < oldest_time) {
            oldest_time = entry->last_used_timestamp;
            oldest_index = i;
        }
    }
    if (oldest_index > -1) {
        entry = _table->entries + oldest_index;
        char qualified[PATH_MAX];
        size_t len = strlen(entry->name);
        assert(len + strlen(_fileprefix) < PATH_MAX);
        strcpy(qualified, _fileprefix);
        strcat(qualified, entry->name);
        printf("Evicting %s\n", entry->name);
        remove(qualified);
        _table->totalbytes -= entry->nbytes;
        if (_table->nentries-- > 1) {
            *entry = _table->entries[_table->nentries];
        }
    }
}

// https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function

static uint64_t _hash(const char* name)
{
    const uint64_t OFFSET = 14695981039346656037ull;
    const uint64_t PRIME = 1099511628211ull;
    const unsigned char* str = (const unsigned char*) name;
    uint64_t hval = OFFSET;
    while (*str) {
        hval *= PRIME;
        hval ^= (uint64_t) *(str++);
    }
    return hval;
}

#undef MIN
#undef MAX
#undef CLAMP
