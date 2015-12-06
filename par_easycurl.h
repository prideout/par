// EASYCURL :: https://github.com/prideout/parg
// Wrapper around libcurl for performing simple synchronous HTTP requests.
//
// The MIT License
// Copyright (c) 2015 Philip Rideout

// -----------------------------------------------------------------------------
// BEGIN PUBLIC API
// -----------------------------------------------------------------------------

typedef unsigned char par_byte;

// Call this before calling any other easycurl function.  The flags are
// currently unused, so you can just pass 0.
void par_easycurl_init(unsigned int flags);

// Allocates a memory buffer and downloads a data blob into it.
// Returns 1 for success and 0 otherwise.  The byte count should be
// pre-allocated.  The caller is responsible for freeing the returned data.
// This does not do any caching!
int par_easycurl_to_memory(char const* url, par_byte** data, int* nbytes);

// Downloads a file from the given URL and saves it to disk.  Returns 1 for
// success and 0 otherwise.
int par_easycurl_to_file(char const* srcurl, char const* dstpath);

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------

#ifdef PAR_EASYCURL_IMPLEMENTATION

#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <curl/curl.h>

static int _ready = 0;

void par_easycurl_init(unsigned int flags)
{
    if (!_ready) {
        curl_global_init(CURL_GLOBAL_SSL);
        _ready = 1;
    }
}

void par_easycurl_shutdown()
{
    if (_ready) {
        curl_global_cleanup();
    }
}

static size_t onheader(void* v, size_t size, size_t nmemb)
{
    size_t n = size * nmemb;
    char* h = v;
    if (n > 14 && !strncasecmp("Last-Modified:", h, 14)) {
        char const* s = h + 14;
        time_t r = curl_getdate(s, 0);
        if (r != -1) {
            // TODO handle last-modified
        }
    } else if (n > 5 && !strncasecmp("ETag:", h, 5)) {
        // TODO handle etag
    }
    return n;
}

typedef struct {
    par_byte* data;
    int nbytes;
} par_easycurl_buffer;

static size_t onwrite(char* contents, size_t size, size_t nmemb, void* udata)
{
    size_t realsize = size * nmemb;
    par_easycurl_buffer* mem = (par_easycurl_buffer*) udata;
    mem->data = realloc(mem->data, mem->nbytes + realsize + 1);
    if (!mem->data) {
        return 0;
    }
    memcpy(mem->data + mem->nbytes, contents, realsize);
    mem->nbytes += realsize;
    mem->data[mem->nbytes] = 0;
    return realsize;
}

#if IOS_EXAMPLE
bool curlToMemory(char const* url, uint8_t** data, int* nbytes)
{
    NSString* nsurl =
    [NSString stringWithCString:url encoding:NSASCIIStringEncoding];
    NSMutableURLRequest* request =
    [NSMutableURLRequest requestWithURL:[NSURL URLWithString:nsurl]];
    [request setTimeoutInterval : TIMEOUT_SECONDS];
    NSURLResponse* response = nil;
    NSError* error = nil;
    // Use the simple non-async API because we're in a secondary thread anyway.
    NSData* nsdata = [NSURLConnection sendSynchronousRequest:request
    returningResponse:&response
    error:&error];
    if (error == nil) {
        *nbytes = (int) [nsdata length];
        *data = (uint8_t*) malloc([nsdata length]);
        memcpy(*data, [nsdata bytes], [nsdata length]);
        return true;
    }
    BLAZE_ERROR("%s\n", [[error localizedDescription] UTF8String]);
    return false;
}

#endif

int par_easycurl_to_memory(char const* url, par_byte** data, int* nbytes)
{
    char errbuf[CURL_ERROR_SIZE] = {0};
    par_easycurl_buffer buffer = {malloc(1), 0};
    long code = 0;
    long status = 0;
    CURL* handle = curl_easy_init();
    curl_easy_setopt(handle, CURLOPT_NOSIGNAL, 1);
    curl_easy_setopt(handle, CURLOPT_ENCODING, "gzip, deflate");
    curl_easy_setopt(handle, CURLOPT_FOLLOWLOCATION, 1);
    curl_easy_setopt(handle, CURLOPT_MAXREDIRS, 8);
    curl_easy_setopt(handle, CURLOPT_FAILONERROR, 1);
    curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, onwrite);
    curl_easy_setopt(handle, CURLOPT_WRITEDATA, &buffer);
    curl_easy_setopt(handle, CURLOPT_HEADERFUNCTION, onheader);
    curl_easy_setopt(handle, CURLOPT_URL, url);
    curl_easy_setopt(handle, CURLOPT_TIMECONDITION, CURL_TIMECOND_IFMODSINCE);
    curl_easy_setopt(handle, CURLOPT_TIMEVALUE, 0);
    curl_easy_setopt(handle, CURLOPT_HTTPHEADER, 0);
    curl_easy_setopt(handle, CURLOPT_TIMEOUT, 60);
    curl_easy_setopt(handle, CURLOPT_ERRORBUFFER, errbuf);
    curl_easy_setopt(handle, CURLOPT_SSL_VERIFYHOST, 0);
    CURLcode res = curl_easy_perform(handle);
    if (res != CURLE_OK) {
        printf("CURL Error: %s\n", errbuf);
        return 0;
    }
    curl_easy_getinfo(handle, CURLINFO_CONDITION_UNMET, &code);
    curl_easy_getinfo(handle, CURLINFO_RESPONSE_CODE, &status);
    if (status == 304 || status >= 400) {
        return 0;
    }
    *data = buffer.data;
    *nbytes = buffer.nbytes;
    curl_easy_cleanup(handle);
    return 1;
}

int par_easycurl_to_file(char const* srcurl, char const* dstpath)
{
    long code = 0;
    long status = 0;
    FILE* filehandle = fopen(dstpath, "wb");
    if (!filehandle) {
        printf("Unable to open %s for writing.\n", dstpath);
        return 0;
    }
    CURL* handle = curl_easy_init();
    curl_easy_setopt(handle, CURLOPT_NOSIGNAL, 1);
    curl_easy_setopt(handle, CURLOPT_ENCODING, "gzip, deflate");
    curl_easy_setopt(handle, CURLOPT_FOLLOWLOCATION, 1);
    curl_easy_setopt(handle, CURLOPT_MAXREDIRS, 8);
    curl_easy_setopt(handle, CURLOPT_FAILONERROR, 1);
    curl_easy_setopt(handle, CURLOPT_WRITEDATA, filehandle);
    curl_easy_setopt(handle, CURLOPT_HEADERFUNCTION, onheader);
    curl_easy_setopt(handle, CURLOPT_URL, srcurl);
    curl_easy_setopt(handle, CURLOPT_TIMECONDITION, CURL_TIMECOND_IFMODSINCE);
    curl_easy_setopt(handle, CURLOPT_TIMEVALUE, 0);
    curl_easy_setopt(handle, CURLOPT_HTTPHEADER, 0);
    curl_easy_setopt(handle, CURLOPT_TIMEOUT, 60);
    curl_easy_perform(handle);
    curl_easy_getinfo(handle, CURLINFO_CONDITION_UNMET, &code);
    curl_easy_getinfo(handle, CURLINFO_RESPONSE_CODE, &status);
    fclose(filehandle);
    if (status == 304 || status >= 400) {
        remove(dstpath);
        return 0;
    }
    curl_easy_cleanup(handle);
    return 1;
}

#endif
