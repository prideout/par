[![Build Status](https://travis-ci.org/prideout/par.svg?branch=master)](https://travis-ci.org/prideout/par)

## par

Single-file C libraries under the MIT license.

library    | description
------------------- | ----
**par_msquares.h** | efficient marching squares implementation
**par_easycurl.h** | simple HTTP requests using libcurl
**par_filecache.h** | LRU caching on your device's filesystem
**par_bluenoise.h** | generate progressive 2D point sequences

## tests

To run tests, you need CMake and libcurl.  On OS X, these can be installed with homebrew:

```bash
$ brew install cmake pkg-config curl
```

Here's how you can tell CMake to use the CMakeLists in the `test` folder, placing all the messy stuff in a new folder called `build`.

```bash
$ cmake test -Bbuild   # Create makefiles
$ cmake --build build  # Invoke the build
```

The tests are executed by simply running the programs:
```bash
$ build/test_msquares
$ build/test_bluenoise
```