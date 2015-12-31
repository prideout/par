#include "describe.h"

#define PAR_SHAPES_IMPLEMENTATION
#include "par_shapes.h"

#include <fcntl.h>
#include <unistd.h>

static int fileexists(const char* filename)
{
    return access(filename, F_OK) != -1;
}

int main()
{
    describe("par_shapes_list_parametric") {
        it("should return several reasonably-sized strings") {
            char const * const * list = par_shapes_list_parametric();
            int nshapes = 0;
            while (*list) {
                const char* name = *list;
                assert(strlen(name) > 0 && strlen(name) < 128);
                nshapes++;
                list++;
            }
            assert(nshapes == 5);
        }
    }

    describe("par_shapes_create_parametric") {
        it("should fail when stacks or slices are less than 3") {
            par_shapes_mesh* bad1 = par_shapes_create_parametric(
                "sphere", 2, 3, 0);
            par_shapes_mesh* bad2 = par_shapes_create_parametric(
                "sphere", 3, 2, 0);
            par_shapes_mesh* good = par_shapes_create_parametric(
                "sphere", 3, 3, 0);
            assert_null(bad1);
            assert_null(bad2);
            assert_ok(good);
            par_shapes_free(good);
        }
        it("should fail with a bogus string") {
            par_shapes_mesh const* bad = par_shapes_create_parametric(
                "bogus", 3, 3, 0);
            assert_null(bad);
        }
        it("should generate correct number of faces and vertices") {
            par_shapes_mesh* m = par_shapes_create_parametric(
                "sphere", 5, 6, 0);
            assert_equal(m->npoints, 42);
            assert_equal(m->ntriangles, 60);
            par_shapes_free(m);
        }
    }

    describe("par_shapes_export") {
        it("should generate an OBJ file") {
            par_shapes_mesh* m;
            m = par_shapes_create_parametric("sphere", 5, 6, 0);
            par_shapes_export(m, "test_shapes_sphere.obj");
            assert_ok(fileexists("test_shapes_sphere.obj"));
            par_shapes_free(m);
            m = par_shapes_create_parametric("plane", 5, 6, 0);
            par_shapes_export(m, "test_shapes_plane.obj");
            assert_ok(fileexists("test_shapes_plane.obj"));
            par_shapes_free(m);
            m = par_shapes_create_parametric("cylinder", 5, 20, 0);
            par_shapes_export(m, "test_shapes_cylinder.obj");
            assert_ok(fileexists("test_shapes_cylinder.obj"));
            par_shapes_free(m);
            m = par_shapes_create_parametric("torus", 7, 10, 0);
            par_shapes_export(m, "test_shapes_torus.obj");
            assert_ok(fileexists("test_shapes_torus.obj"));
            par_shapes_free(m);
            m = par_shapes_create_parametric("klein", 10, 20, 0);
            par_shapes_export(m, "test_shapes_klein.obj");
            assert_ok(fileexists("test_shapes_klein.obj"));
            par_shapes_free(m);
        }
    }

    describe("par_shapes_merge") {
        it("concatenate two meshes") {
            par_shapes_mesh* a, *b;
            a = par_shapes_create_parametric("klein", 10, 20, 0);
            int npts = a->npoints;
            int ntris = a->ntriangles;
            b = par_shapes_create_parametric("plane", 3, 3, 0);
            par_shapes_merge(a, b);
            assert_equal(a->npoints, npts + b->npoints);
            assert_equal(a->ntriangles, ntris + b->ntriangles);
            par_shapes_export(a, "test_shapes_merged.obj");
            par_shapes_free(a);
            par_shapes_free(b);
        }
    }

    return assert_failures();
}
