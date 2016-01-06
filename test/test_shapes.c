#include "describe.h"

#define PAR_SHAPES_IMPLEMENTATION
#include "par_shapes.h"

#include <fcntl.h>
#include <unistd.h>

#define PAR_PI (3.14159265359)

int main()
{
    describe("par_shapes_export") {
        it("should generate an OBJ file") {
            par_shapes_mesh* m;

            m = par_shapes_create_cylinder(50, 20);
            par_shapes_export(m, "build/test_shapes_cylinder.obj");
            par_shapes_free(m);

            m = par_shapes_create_torus(7, 10, 0.1);
            par_shapes_export(m, "build/test_shapes_torus.obj");
            par_shapes_free(m);

            m = par_shapes_create_parametric_sphere(5, 6);
            par_shapes_mesh* m2 = par_shapes_weld(m, 0.001, NULL);
            par_shapes_export(m2, "build/test_shapes_psphere.obj");
            par_shapes_free(m);
            par_shapes_free(m2);

            m = par_shapes_create_subdivided_sphere(2);
            par_shapes_export(m, "build/test_shapes_ssphere.obj");
            par_shapes_free(m);

            m = par_shapes_create_klein_bottle(10, 20);
            par_shapes_export(m, "build/test_shapes_klein.obj");
            par_shapes_free(m);

            m = par_shapes_create_trefoil_knot(20, 100, 0.2);
            par_shapes_export(m, "build/test_shapes_trefoil.obj");
            par_shapes_free(m);

            m = par_shapes_create_hemisphere(5, 6);
            par_shapes_export(m, "build/test_shapes_hemisphere.obj");
            par_shapes_free(m);

            m = par_shapes_create_plane(5, 6);
            par_shapes_export(m, "build/test_shapes_plane.obj");
            par_shapes_free(m);

            m = par_shapes_create_icosahedron();
            par_shapes_export(m, "build/test_shapes_icosahedron.obj");
            par_shapes_free(m);

            m = par_shapes_create_dodecahedron();
            par_shapes_export(m, "build/test_shapes_dodecahedron.obj");
            par_shapes_free(m);

            m = par_shapes_create_octohedron();
            par_shapes_export(m, "build/test_shapes_octohedron.obj");
            par_shapes_free(m);

            m = par_shapes_create_tetrahedron();
            par_shapes_export(m, "build/test_shapes_tetrahedron.obj");
            par_shapes_free(m);

            m = par_shapes_create_cube();
            par_shapes_export(m, "build/test_shapes_cube.obj");
            par_shapes_free(m);

            m = par_shapes_create_rock(1, 3);
            par_shapes_export(m, "build/test_shapes_rock.obj");
            par_shapes_free(m);

            float center[3] = {0, 0, 0};
            float normal[3] = {0, 0, 1};
            m = par_shapes_create_disk(1, 5, center, normal);
            par_shapes_export(m, "build/test_shapes_disk.obj");
            par_shapes_free(m);

        }
    }

    describe("par_shapes_create_cylinder") {
        it("should fail when the number of stacks or slices is invalid") {
            par_shapes_mesh* bad1 = par_shapes_create_cylinder(1, 1);
            par_shapes_mesh* bad2 = par_shapes_create_cylinder(1, 3);
            par_shapes_mesh* good = par_shapes_create_cylinder(3, 1);
            assert_null(bad1);
            assert_null(bad2);
            assert_ok(good);
            par_shapes_free(good);
        }
        it("should generate correct number of faces and vertices") {
            par_shapes_mesh* m = par_shapes_create_cylinder(5, 6);
            assert_equal(m->npoints, 42);
            assert_equal(m->ntriangles, 60);
            par_shapes_free(m);
        }
    }

    describe("par_shapes_merge") {
        it("should concatenate two meshes") {
            par_shapes_mesh* a, *b;
            a = par_shapes_create_klein_bottle(10, 20);
            int npts = a->npoints;
            int ntris = a->ntriangles;
            b = par_shapes_create_plane(3, 3);
            par_shapes_merge(a, b);
            assert_equal(a->npoints, npts + b->npoints);
            assert_equal(a->ntriangles, ntris + b->ntriangles);
            par_shapes_free(a);
            par_shapes_free(b);
        }
    }

    describe("transforms") {
        it("should support translation") {
            par_shapes_mesh* a, *b;
            a = par_shapes_create_cylinder(20, 3);
            b = par_shapes_create_cylinder(4, 3);
            par_shapes_translate(a, 0.5, 0.5, 0.25);
            par_shapes_merge(a, b);
            par_shapes_free(a);
            par_shapes_free(b);
        }
        it("should support rotation") {
            par_shapes_mesh* a, *b;
            a = par_shapes_create_cylinder(20, 3);
            b = par_shapes_create_cylinder(4, 3);
            float axis1[3] = {0, 1, 0};
            float axis2[3] = {0, 0, 1};
            par_shapes_rotate(a, PAR_PI * 0.5, axis1);
            par_shapes_rotate(a, PAR_PI * 0.25, axis2);
            par_shapes_merge(a, b);
            par_shapes_free(a);
            par_shapes_free(b);
        }
        it("should support non-uniform scale") {
            par_shapes_mesh* a;
            a = par_shapes_create_cylinder(15, 3);
            par_shapes_scale(a, 1, 1, 5);
            par_shapes_free(a);
        }
    }

    describe("misc shapes") {
        it("create an orientable disk in 3-space") {
            int slices = 32;
            float aradius = 1;
            float anormal[3] = {0, 0, 1};
            float acenter[3] = {0, 0, 0};
            par_shapes_mesh* a, *b;
            a = par_shapes_create_disk(aradius, slices, acenter, anormal);
            float bradius = 0.2;
            float bcenter[3] = {0, 0, 0.2};
            float bnormal[3] = {0, 1, 0};
            b = par_shapes_create_disk(bradius, slices, bcenter, bnormal);
            par_shapes_merge(a, b);
            par_shapes_free(a);
            par_shapes_free(b);
        }
        it("create a rock on the Y plane") {
            int slices = 32;
            float radius = 2;
            float normal[3] = {0, 1, 0};
            float center[3] = {0, 0, 0};
            par_shapes_mesh* a, *b;
            a = par_shapes_create_disk(radius, slices, center, normal);
            b = par_shapes_create_rock(1, 2);
            float aabb[6];
            par_shapes_compute_aabb(b, aabb);
            par_shapes_translate(b, 0, -aabb[1] / 2, 0);
            par_shapes_merge(a, b);
            par_shapes_free(a);
            par_shapes_free(b);
        }
        it("create a polyhedron on the Y plane") {
            int slices = 32;
            float radius = 2;
            float normal[3] = {0, 1, 0};
            float center[3] = {0, 0, 0};
            par_shapes_mesh* a, *b;
            a = par_shapes_create_disk(radius, slices, center, normal);
            b = par_shapes_create_dodecahedron();
            par_shapes_translate(b, 0, 0.934, 0);
            par_shapes_merge(a, b);
            par_shapes_free(a);
            par_shapes_free(b);
        }
        it("create a rounded cylinder via composition") {
            const float O[3] = {0, 0, 0};
            const float I[3] = {1, 0, 0};
            const float J[3] = {0, 1, 0};
            const float K[3] = {0, 0, 1};
            const float top_center[3] = {0, 1.2, 0};
            const int tess = 30;
            par_shapes_mesh *a, *b, *c, *d;
            a = par_shapes_create_disk(2.5, tess, O, J);
            b = par_shapes_create_cylinder(tess, 3);
            c = par_shapes_create_torus(15, tess, 0.1);
            d = par_shapes_create_disk(1, tess, top_center, J);
            par_shapes_rotate(c, PAR_PI / tess, K);
            par_shapes_translate(c, 0, 0, 1);
            par_shapes_scale(b, 1.2, 1.2, 1);
            par_shapes_merge(b, c);
            par_shapes_rotate(b, -PAR_PI * 0.5, I);
            par_shapes_merge(b, d);
            par_shapes_merge(b, a);
            par_shapes_scale(b, 1, 2, 1);
            par_shapes_free(a);
            par_shapes_free(b);
            par_shapes_free(c);
            par_shapes_free(d);
        }
    }

    return assert_failures();
}
