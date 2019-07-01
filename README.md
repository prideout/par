[![Build Status](https://travis-ci.org/prideout/par.svg?branch=master)](https://travis-ci.org/prideout/par)

## par

Single-file C libraries under the MIT license.  Documentation can be found at the top of each header file, but some libraries have an accompanying blog post.

library    | description  | link
------------------- | ---- | ---
**par_streamlines.h** | triangulate wide lines and curves | [blog post](https://prideout.net/blog/par_streamlines/)
**par_shapes.h** | generate parametric surfaces and other simple shapes | [blog post](https://prideout.net/shapes)
**par_shaders.h** | string extraction and concatenation |
**par_sprune.h** | efficient broad-phase collision detection in 2D | [web demo](https://prideout.net/d3cpp/)
**par_easings.h** | Robert Penner's easing functions |
**par_bubbles.h** | pack circles into hierarchical diagrams | [blog post](https://prideout.net/bubbles)
**par_msquares.h** | efficient marching squares implementation | [blog post](https://prideout.net/marching-squares)
**par_bluenoise.h** | generate progressive 2D point sequences | [blog post](https://prideout.net/recursive-wang-tiles)
**par_easycurl.h** | simple HTTP requests using libcurl |
**par_filecache.h** | LRU caching on your device's filesystem |

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
$ build/test_bubbles
$ build/test_shapes
```

## code formatting

This library's code style is strictly enforced to be vertically dense (no consecutive newlines) and horizontally narrow (80 columns or less).

The `tools/format.py` script invokes a two-step code formatting process:

1. Runs `uncrustify` with our custom configuration.  This auto-formats all code in the root folder, up to a point.  For example, it does not enforce the 80-character line constraint because line breaking is best done by a human.
1. Checks for violations that are not otherwise enforced with uncrustify.

The aforementioned Python script is also invoked from Travis, but using the `--check` option, which checks for conformance without editing the code.

Beyond what our uncrustify configuration enforces, the Python script does the following:

- Checks that no lines are more than 80 chars.
- Checks for extra newlines before an end brace.
