extern "C" {
#include "describe.h"
}

#define PAR_OCTASPHERE_IMPLEMENTATION
#include "par_octasphere.h"

#define CGLTF_IMPLEMENTATION
#define CGLTF_WRITE_IMPLEMENTATION
#include "cgltf_write.h"

static void generate_gltf(const char* gltf_path, const char* bin_path, int num_vertices,
        int num_indices, float minpos[3], float maxpos[3]) {
    cgltf_image images[3] = {
        {(char*)"octasphere color", (char*)"octasphere_color.png"},
        {(char*)"octasphere orm", (char*)"octasphere_orm.png"},
        {(char*)"octasphere normal", (char*)"octasphere_normal.png"},
    };
    cgltf_image& color_image = images[0];
    cgltf_image& orm_image = images[1];
    cgltf_image& normal_image = images[2];

    cgltf_texture textures[3] = {
        {(char*)"octasphere color", &color_image},
        {(char*)"octasphere orm", &orm_image},
        {(char*)"octasphere normal", &normal_image},
    };
    cgltf_texture& color_texture = textures[0];
    cgltf_texture& orm_texture = textures[1];
    cgltf_texture& normal_texture = textures[2];

    cgltf_material materials[1] = {};
    cgltf_material& octasphere_material = materials[0];

    octasphere_material.name = (char*)"tile top";
    octasphere_material.alpha_cutoff = 0.5f;
    octasphere_material.has_pbr_metallic_roughness = true;
    octasphere_material.pbr_metallic_roughness = {
        .base_color_texture = {&color_texture, 0, 1.0f},
        .metallic_roughness_texture = {&orm_texture, 0, 1.0f},
        .base_color_factor = {1, 1, 1, 1},
        .metallic_factor = 1,
        .roughness_factor = 1,
    };
    octasphere_material.occlusion_texture = {&orm_texture, 0, 1.0f};
    octasphere_material.normal_texture = {&normal_texture, 0, 1.0f};

    // GEOMETRY

    cgltf_buffer buffers[1] = {};
    cgltf_mesh meshes[1] = {};
    cgltf_buffer_view buffer_views[4] = {};
    cgltf_accessor accessors[4] = {};

    cgltf_attribute attributes[3] = {};
    cgltf_primitive prims[2] = {};

    cgltf_buffer& octasphere_buffer = buffers[0];

    octasphere_buffer.uri = (char*) bin_path;
    octasphere_buffer.size = num_vertices * 32 + num_indices * 2;

    buffer_views[0].buffer = &octasphere_buffer;
    buffer_views[0].size = num_vertices * sizeof(float) * 3;
    buffer_views[0].type = cgltf_buffer_view_type_vertices;

    buffer_views[1].buffer = &octasphere_buffer;
    buffer_views[1].size = num_vertices * sizeof(float) * 3;
    buffer_views[1].offset = buffer_views[0].offset;
    buffer_views[1].offset += buffer_views[0].size;
    buffer_views[1].type = cgltf_buffer_view_type_vertices;

    buffer_views[2].buffer = &octasphere_buffer;
    buffer_views[2].size = num_vertices * sizeof(float) * 2;
    buffer_views[2].offset = buffer_views[1].offset;
    buffer_views[2].offset += buffer_views[1].size;
    buffer_views[2].type = cgltf_buffer_view_type_vertices;

    buffer_views[3].buffer = &octasphere_buffer;
    buffer_views[3].size = num_indices * sizeof(uint16_t);
    buffer_views[3].offset = buffer_views[2].offset;
    buffer_views[3].offset += buffer_views[2].size;
    buffer_views[3].type = cgltf_buffer_view_type_indices;

    accessors[0].buffer_view = &buffer_views[0];
    accessors[0].component_type = cgltf_component_type_r_32f;
    accessors[0].type = cgltf_type_vec3;
    accessors[0].count = num_vertices;
    accessors[0].has_min = accessors[0].has_max = true;

    for (int i = 0; i < 3; i++) {
        accessors[0].min[i] = minpos[i];
        accessors[0].max[i] = maxpos[i];
    }

    accessors[1].buffer_view = &buffer_views[1];
    accessors[1].component_type = cgltf_component_type_r_32f;
    accessors[1].type = cgltf_type_vec3;
    accessors[1].count = num_vertices;

    accessors[2].buffer_view = &buffer_views[2];
    accessors[2].component_type = cgltf_component_type_r_32f;
    accessors[2].type = cgltf_type_vec2;
    accessors[2].count = num_vertices;

    accessors[3].buffer_view = &buffer_views[3];
    accessors[3].component_type = cgltf_component_type_r_16u;
    accessors[3].type = cgltf_type_scalar;
    accessors[3].count = num_indices;

    // ATTRIBUTES

    attributes[0].name = (char*)"POSITION";
    attributes[0].data = &accessors[0];
    attributes[0].type = cgltf_attribute_type_position;

    attributes[1].name = (char*)"NORMAL";
    attributes[1].data = &accessors[1];
    attributes[1].type = cgltf_attribute_type_normal;

    attributes[2].name = (char*)"TEXCOORD_0";
    attributes[2].data = &accessors[2];
    attributes[2].type = cgltf_attribute_type_texcoord;

    // PRIMITIVES AND MESHES

    cgltf_primitive& prim = prims[0];
    prim.type = cgltf_primitive_type_triangles;
    prim.attributes = attributes;
    prim.attributes_count = 3;
    prim.indices = &accessors[3];
    prim.material = &octasphere_material;

    cgltf_mesh& mesh = meshes[0];
    mesh.name = (char*)"octasphere mesh";
    mesh.primitives = &prim;
    mesh.primitives_count = 1;

    // NODE HIERARCHY

    cgltf_node nodes[2];
    cgltf_node* pnodes[2] = {&nodes[0], &nodes[1]};

    nodes[0] = {
        .name = (char*)"root",
        .parent = {},
        .children = &pnodes[1],
        .children_count = 1,
        .skin = {},
        .mesh = {},
        .camera = {},
        .light = {},
        .weights = {},
        .weights_count = 0,
        .has_translation = false,
        .has_rotation = false,
        .has_scale = false,
        .has_matrix = true,
        .translation = {},
        .rotation = {},
        .scale = {},
        .matrix = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
        },
    };

    nodes[1] = {
        .name = (char*) "octasphere",
        .parent = {},
        .children = {},
        .children_count = 0,
        .skin = {},
        .mesh = &meshes[0]
    };

    cgltf_scene scene = {};
    scene.nodes = pnodes;
    scene.nodes_count = sizeof(pnodes) / sizeof(pnodes[0]);

    // ASSET

    cgltf_data asset = {};
    asset.file_type = cgltf_file_type_gltf;
    asset.asset.generator = (char*)"test_octasphere";
    asset.asset.version = (char*)"2.0";

    asset.buffers = buffers;
    asset.buffers_count = sizeof(buffers) / sizeof(buffers[0]);

    asset.buffer_views = buffer_views;
    asset.buffer_views_count = sizeof(buffer_views) / sizeof(buffer_views[0]);

    asset.accessors = accessors;
    asset.accessors_count = sizeof(accessors) / sizeof(accessors[0]);

    asset.meshes = meshes;
    asset.meshes_count = sizeof(meshes) / sizeof(meshes[0]);

    asset.materials = materials;
    asset.materials_count = sizeof(materials) / sizeof(materials[0]);

    asset.images = images;
    asset.images_count = sizeof(images) / sizeof(images[0]);

    asset.textures = textures;
    asset.textures_count = sizeof(textures) / sizeof(textures[0]);

    asset.nodes = nodes;
    asset.nodes_count = sizeof(nodes) / sizeof(nodes[0]);

    asset.scenes = &scene;
    asset.scenes_count = 1;

    cgltf_options options = {};
    cgltf_write_file(&options, gltf_path, &asset);
}

int main()
{
    describe("octasphere generator") {

        it("should generate a tile-like shape") {
            par_octasphere_config config = {
                .corner_radius = 0.1,
                .width = 1.2,
                .height = 1.2,
                .depth = 0.3,
                .num_subdivisions = 2,
                .uv_mode = PAR_OCTASPHERE_UV_LATLONG,
            };

            uint32_t indices_per_tile;
            uint32_t vertices_per_tile;
            par_octasphere_get_counts(&config, &indices_per_tile, &vertices_per_tile);

            par_octasphere_mesh octasphere = {};
            octasphere.positions = (float*)malloc(vertices_per_tile * 12);
            octasphere.normals = (float*)malloc(vertices_per_tile * 12);
            octasphere.texcoords = (float*)malloc(vertices_per_tile * 8);
            octasphere.indices = (uint16_t*)malloc(indices_per_tile * 2);
            par_octasphere_populate(&config, &octasphere);

            free(octasphere.positions);
            free(octasphere.normals);
            free(octasphere.texcoords);
            free(octasphere.indices);
        }

        it("should generate rounded cube into glTF + bin files") {
            par_octasphere_config config = {
                .corner_radius = 0.4,
                .width = 1.2,
                .height = 1.2,
                .depth = 1.2,
                .num_subdivisions = 3,
                .uv_mode = PAR_OCTASPHERE_UV_LATLONG,
            };

            uint32_t num_indices;
            uint32_t num_vertices;
            par_octasphere_get_counts(&config, &num_indices, &num_vertices);

            par_octasphere_mesh octasphere = {};
            octasphere.positions = (float*)malloc(num_vertices * 12);
            octasphere.normals = (float*)malloc(num_vertices * 12);
            octasphere.texcoords = (float*)malloc(num_vertices * 8);
            octasphere.indices = (uint16_t*)malloc(num_indices * 2);
            par_octasphere_populate(&config, &octasphere);

            float minpos[3] = { 99,  99,  99};
            float maxpos[3] = {-99, -99, -99};
            // TODO: compute minpos and maxpos

            FILE* file = fopen("octasphere.bin", "wb");
            if (!file) {
                puts("unable to open bin file for writing");
                exit(1);
            }
            fwrite(octasphere.positions, 12, num_vertices, file);
            fwrite(octasphere.normals, 12, num_vertices, file);
            fwrite(octasphere.texcoords, 8, num_vertices, file);
            fwrite(octasphere.indices, 2, num_indices, file);
            fclose(file);

            free(octasphere.positions);
            free(octasphere.normals);
            free(octasphere.texcoords);
            free(octasphere.indices);

            generate_gltf("octasphere.gltf", "octasphere.bin", num_vertices, num_indices,
                    minpos, maxpos);
        }
    }

    return assert_failures();
}
