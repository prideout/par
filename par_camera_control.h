// CAMERA CONTROL :: https://github.com/prideout/par
// Enables orbit controls (a.k.a. tumble, arcball, trackball) or pan-and-zoom like Google Maps.
//
// This simple library controls a camera that orbits or pans over a 3D object or terrain. No
// assumptions are made about the renderer or platform. In a sense, this is just a math library.
// Clients notify the controller of generic input events (e.g. grab_begin, grab_move, grab_end)
// and retrieve the look-at vectors (position, target, up) or 4x4 matrices for the camera.
//
// In map mode, users can control their viewing position by grabbing and dragging locations in the
// scene. Sometimes this is known as "through-the-lens" camera control. In this mode the controller
// takes an optional raycast callback to support precise grabbing behavior. If this is not required
// for your use case (e.g. a top-down terrain with an orthgraphic projection), provide NULL for the
// callback and the library will simply raycast against the ground plane.
//
// When the controller is in orbit mode, the orientation of the camera is defined by a Y-axis
// rotation followed by an X-axis rotation. Additionally, the camera can fly forward or backward
// along the viewing direction.
//
// For a complex usage example, go to:
// https://github.com/prideout/camera_demo
//
// Distributed under the MIT License, see bottom of file.

#ifndef PAR_CAMERA_CONTROL_H
#define PAR_CAMERA_CONTROL_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PARCC_USE_DOUBLE
typedef double parcc_float;
#else
typedef float parcc_float;
#endif

// Opaque handle to a camera controller.
typedef struct parcc_context_s parcc_context;

// The camera controller can be configured using either a VERTICAL or HORIZONTAL field of view.
// This specifies which of the two FOV angles should be held constant. For example, if you use a
// horizontal FOV, shrinking the viewport width will change the height of the frustum, but will
// leave the frustum width intact.
typedef enum {
    PARCC_VERTICAL,
    PARCC_HORIZONTAL,
} parcc_fov;

// The controller can be configured in orbit mode or pan-and-zoom mode.
typedef enum {
    PARCC_ORBIT,  // aka tumble, trackball, or arcball
    PARCC_MAP,    // pan and zoom like Google Maps
} parcc_mode;

// Pan and zoom constraints for MAP mode.
typedef enum {
    // No constraints except that map_min_distance is enforced.
    PARCC_CONSTRAIN_NONE,

    // Constrains pan and zoom to limit the viewport's extent along the FOV axis so that it always
    // lies within the map_extent. With this constraint, it is possible to see the entire map at
    // once, but some portion of the map must always be visible.
    PARCC_CONSTRAIN_AXIS,

    // Constrains pan and zoom to limit the viewport's extent into the map_extent. With this
    // constraint, it may be impossible to see the entire map at once, but users can never see any
    // of the empty void that lies outside the map extent.
    PARCC_CONSTRAIN_FULL,
} parcc_constraint;

// Optional user-provided ray casting function to enable precise panning behavior.
typedef bool (*parcc_raycast_fn)(const parcc_float origin[3], const parcc_float dir[3],
    parcc_float* t, void* userdata);

// The parcc_properties structure holds all user-controlled state in the library.
// Many fields are swapped with fallback values values if they are zero-filled.

typedef struct {
    // REQUIRED PROPERTIES
    parcc_mode mode;         // must be PARCC_ORBIT or PARCC_MAP
    int viewport_width;      // horizontal extent in pixels
    int viewport_height;     // vertical extent in pixels
    parcc_float near_plane;  // distance between camera and near clipping plane
    parcc_float far_plane;   // distance between camera and far clipping plane

    // PROPERTIES WITH DEFAULT VALUES
    parcc_fov fov_orientation;   // defaults to PARCC_VERTICAL
    parcc_float fov_degrees;     // full field-of-view angle (not half-angle), defaults to 33.
    parcc_float zoom_speed;      // defaults to 0.01
    parcc_float home_target[3];  // world-space coordinate, defaults to (0,0,0)
    parcc_float home_upward[3];  // unit-length vector, defaults to (0,1,0)

    // MAP-MODE PROPERTIES
    parcc_float map_extent[2];          // (required) size of quad centered at home_target
    parcc_float map_plane[4];           // plane equation with normalized XYZ, defaults to (0,0,1,0)
    parcc_constraint map_constraint;    // defaults to PARCC_CONSTRAIN_NONE
    parcc_float map_min_distance;       // constrains zoom using distance between camera and plane
    parcc_raycast_fn raycast_function;  // defaults to a simple plane intersector
    void* raycast_userdata;             // arbitrary data for the raycast callback

    // ORBIT-MODE PROPERTIES
    parcc_float home_vector[3];    // (required) vector from home_target to initial eye position
    parcc_float orbit_speed[2];    // rotational speed (defaults to 0.01)
    parcc_float orbit_zoom_speed;  // zoom speed (defaults to 0.01)

} parcc_properties;

// The parcc_frame structure holds captured camera state for Van Wijk animation and bookmarks.
// From the user's perspective, this should be treated as an opaque structure.
// clang-format off

typedef struct {
    parcc_mode mode;
    union {
        struct { parcc_float extent, center[2]; };
        struct { parcc_float phi, theta, pivot_distance, pivot[3]; };
    };
} parcc_frame;

// clang-format on
// CONTROLLER CONSTRUCTOR AND DESTRUCTOR

// The constructor is the only function in the library that performs heap allocation. It
// does not retain the given properties pointer, it simply copies values out of it.

parcc_context* parcc_create_context(const parcc_properties* props);
void parcc_destroy_context(parcc_context* ctx);

// PROPERTY SETTERS AND GETTERS

// The client owns its own instance of the property struct and these functions simply copy values in
// or out of the given struct. Changing some properties might cause a small amount of work to be
// performed.

void parcc_set_properties(parcc_context* context, const parcc_properties* props);
void parcc_get_properties(const parcc_context* context, parcc_properties* out);

// CAMERA RETRIEVAL FUNCTIONS

void parcc_get_look_at(const parcc_context* ctx, parcc_float eyepos[3], parcc_float target[3],
    parcc_float upward[3]);
void parcc_get_matrices(const parcc_context* ctx, parcc_float projection[16], parcc_float view[16]);

// SCREEN-SPACE FUNCTIONS FOR USER INTERACTION

// Each of these functions take winx / winy coords.
// - The winx coord should be in [0, viewport_width) where 0 is the left-most column.
// - The winy coord should be in [0, viewport_height) where 0 is the top-most row.
//
// The scrolldelta argument is used for zooming. Positive values indicate "zoom in" in MAP mode or
// "move forward" in ORBIT mode. This gets scaled by zoom_speed. In MAP mode, the zoom speed is also
// scaled by distance-to-ground. To prevent zooming in too far, use a non-zero value for
// map_min_distance.
//
// The strafe argument exists only for ORBIT mode and is typically associated with the right mouse
// button or two-finger dragging. This is used to pan the view. Note that orbit mode maintains a
// "pivot point" which is initially set to home_target. When flying forward or backward, the pivot
// does not move. However strafing will cause it to move around. This matches sketchfab behavior.
// When flying past the orbit point, the controller enters a "flipped" state to prevent the flight
// direction from suddenly changing.

void parcc_grab_begin(parcc_context* context, int winx, int winy, bool strafe);
void parcc_grab_update(parcc_context* context, int winx, int winy);
void parcc_grab_end(parcc_context* context);
void parcc_zoom(parcc_context* context, int winx, int winy, parcc_float scrolldelta);
bool parcc_raycast(parcc_context* context, int winx, int winy, parcc_float result[3]);

// BOOKMARKING AND VAN WIJK INTERPOLATION FUNCTIONS

parcc_frame parcc_get_current_frame(const parcc_context* context);
parcc_frame parcc_get_home_frame(const parcc_context* context);
void parcc_goto_frame(parcc_context* context, parcc_frame state);
parcc_frame parcc_interpolate_frames(parcc_frame a, parcc_frame b, double t);
double parcc_get_interpolation_duration(parcc_frame a, parcc_frame b);

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// END PUBLIC API
// -----------------------------------------------------------------------------
#ifdef PAR_CAMERA_CONTROL_IMPLEMENTATION

#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>

#define PARCC_PI (3.14159265359)
#define PARCC_MIN(a, b) (a > b ? b : a)
#define PARCC_MAX(a, b) (a > b ? a : b)
#define PARCC_CLAMP(v, lo, hi) PARCC_MAX(lo, PARCC_MIN(hi, v))
#define PARCC_CALLOC(T, N) ((T*)calloc(N * sizeof(T), 1))
#define PARCC_FREE(BUF) free(BUF)
#define PARCC_SWAP(T, A, B) \
    {                       \
        T tmp = B;          \
        B = A;              \
        A = tmp;            \
    }

inline void parcc_float4_set(parcc_float dst[4], parcc_float x, parcc_float y, parcc_float z,
    parcc_float w) {
    dst[0] = x;
    dst[1] = y;
    dst[2] = z;
    dst[3] = w;
}

inline parcc_float parcc_float4_dot(const parcc_float a[4], const parcc_float b[4]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

inline void parcc_float3_set(parcc_float dst[3], parcc_float x, parcc_float y, parcc_float z) {
    dst[0] = x;
    dst[1] = y;
    dst[2] = z;
}

inline void parcc_float3_add(parcc_float dst[3], const parcc_float a[3], const parcc_float b[3]) {
    dst[0] = a[0] + b[0];
    dst[1] = a[1] + b[1];
    dst[2] = a[2] + b[2];
}

inline void parcc_float3_macc(parcc_float dst[3], const parcc_float src[3], parcc_float scale) {
    dst[0] += src[0] * scale;
    dst[1] += src[1] * scale;
    dst[2] += src[2] * scale;
}

inline void parcc_float3_subtract(parcc_float dst[3], const parcc_float a[3],
    const parcc_float b[3]) {
    dst[0] = a[0] - b[0];
    dst[1] = a[1] - b[1];
    dst[2] = a[2] - b[2];
}

inline parcc_float parcc_float3_dot(const parcc_float a[3], const parcc_float b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void parcc_float3_cross(parcc_float dst[3], const parcc_float a[3], const parcc_float b[3]) {
    dst[0] = a[1] * b[2] - a[2] * b[1];
    dst[1] = a[2] * b[0] - a[0] * b[2];
    dst[2] = a[0] * b[1] - a[1] * b[0];
}

inline void parcc_float3_scale(parcc_float dst[3], parcc_float v) {
    dst[0] *= v;
    dst[1] *= v;
    dst[2] *= v;
}

inline parcc_float parcc_float3_length(const parcc_float dst[3]) {
    return sqrtf(parcc_float3_dot(dst, dst));
}

inline void parcc_float3_normalize(parcc_float dst[3]) {
    parcc_float3_scale(dst, 1.0f / parcc_float3_length(dst));
}

inline void parcc_float3_copy(parcc_float dst[3], const parcc_float src[3]) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
}

inline void parcc_float3_copy_to_vec4(parcc_float dst[4], const parcc_float src[3], parcc_float w) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = w;
}

inline void parcc_float3_lerp(parcc_float dst[3], const parcc_float a[3], const parcc_float b[3],
    parcc_float t) {
    dst[0] = a[0] * (1 - t) + b[0] * t;
    dst[1] = a[1] * (1 - t) + b[1] * t;
    dst[2] = a[2] * (1 - t) + b[2] * t;
}

inline void parcc_float16_look_at(float dst[16], const float eye[3], const float target[3],
    const float up[3]) {
    parcc_float v3X[3];
    parcc_float v3Y[3];
    parcc_float v3Z[3];

    parcc_float3_copy(v3Y, up);
    parcc_float3_normalize(v3Y);

    parcc_float3_subtract(v3Z, eye, target);
    parcc_float3_normalize(v3Z);

    parcc_float3_cross(v3X, v3Y, v3Z);
    parcc_float3_normalize(v3X);

    parcc_float3_cross(v3Y, v3Z, v3X);

    parcc_float4_set(dst + 0, v3X[0], v3Y[0], v3Z[0], 0);
    parcc_float4_set(dst + 4, v3X[1], v3Y[1], v3Z[1], 0);
    parcc_float4_set(dst + 8, v3X[2], v3Y[2], v3Z[2], 0);
    parcc_float4_set(dst + 12,                     //
        -parcc_float3_dot(v3X, eye),               //
        -parcc_float3_dot(v3Y, eye),               //
        -parcc_float3_dot(v3Z, eye), 1.0);
}

inline void parcc_float16_perspective_y(float dst[16], float fovy_degrees, float aspect_ratio,
    float near, float far) {
    const parcc_float fovy_radians = fovy_degrees * PARCC_PI / 180;
    const parcc_float f = tan(PARCC_PI / 2.0 - 0.5 * fovy_radians);
    const parcc_float rangeinv = 1.0f / (near - far);
    dst[0] = f / aspect_ratio;
    dst[1] = 0;
    dst[2] = 0;
    dst[3] = 0;
    dst[4] = 0;
    dst[5] = f;
    dst[6] = 0;
    dst[7] = 0;
    dst[8] = 0;
    dst[9] = 0;
    dst[10] = (near + far) * rangeinv;
    dst[11] = -1;
    dst[12] = 0;
    dst[13] = 0;
    dst[14] = ((near * far) * rangeinv) * 2.0f;
    dst[15] = 0;
}

inline void parcc_float16_perspective_x(float dst[16], float fovy_degrees, float aspect_ratio,
    float near, float far) {
    const parcc_float fovy_radians = fovy_degrees * PARCC_PI / 180;
    const parcc_float f = tan(PARCC_PI / 2.0 - 0.5 * fovy_radians);
    const parcc_float rangeinv = 1.0 / (near - far);
    dst[0] = f;
    dst[1] = 0;
    dst[2] = 0;
    dst[3] = 0;
    dst[4] = 0;
    dst[5] = f * aspect_ratio;
    dst[6] = 0;
    dst[7] = 0;
    dst[8] = 0;
    dst[9] = 0;
    dst[10] = (near + far) * rangeinv;
    dst[11] = -1;
    dst[12] = 0;
    dst[13] = 0;
    dst[14] = ((near * far) * rangeinv) * 2.0;
    dst[15] = 0;
}

// Implementation note about the "parcc_frame" POD. This is an abbreviated camera state
// used for animation and bookmarking.
//
// MAP mode:
// - zoom level is represented with the extent of the rectangle formed by the intersection of
//   the frustum with the viewing plane at home_target. It is either a width or a height, depending
//   on fov_orientation.
// - the pan offset is stored as a 2D vector from home_target that gets projected to map_plane.
//
// ORBIT mode:
// - phi   = X-axis rotation in [-pi/2, +pi/2] (applies first)
// - theta = Y-axis rotation in [-pi, +pi]     (applies second)
// - pivot is initialized to home_center but might be changed via strafe
// - pivot_distance is the distance between eye and pivot (negative distance = orbit_flipped)

typedef enum { PARCC_GRAB_NONE, PARCC_GRAB, PARCC_GRAB_STRAFE } parcc_grab_state;

static const parcc_float PARCC_MAX_PHI = PARCC_PI / 2.0 - 0.001;

struct parcc_context_s {
    parcc_properties props;
    parcc_float eyepos[3];
    parcc_float target[3];
    parcc_grab_state grabbing;
    parcc_float grab_point_pivot[3];
    parcc_float grab_point_far[3];
    parcc_float grab_point_world[3];
    parcc_float grab_point_eyepos[3];
    parcc_float grab_point_target[3];
    parcc_frame grab_frame;
    int grab_winx;
    int grab_winy;
    parcc_float orbit_pivot[3];
    bool orbit_flipped;
};

static bool parcc_raycast_plane(const parcc_float origin[3], const parcc_float dir[3],
    parcc_float* t, void* userdata);

static void parcc_get_ray_far(parcc_context* context, int winx, int winy, parcc_float result[3]);

static void parcc_move_with_constraints(parcc_context* context, const parcc_float eyepos[3],
    const parcc_float target[3]);

parcc_context* parcc_create_context(const parcc_properties* props) {
    parcc_context* context = PARCC_CALLOC(parcc_context, 1);
    parcc_set_properties(context, props);
    parcc_goto_frame(context, parcc_get_home_frame(context));
    return context;
}

void parcc_get_properties(const parcc_context* context, parcc_properties* props) {
    *props = context->props;
}

void parcc_set_properties(parcc_context* context, const parcc_properties* pprops) {
    parcc_properties props = *pprops;
    if (props.fov_degrees == 0) {
        props.fov_degrees = 33;
    }
    if (props.zoom_speed == 0) {
        props.zoom_speed = 0.01;
    }
    if (parcc_float3_dot(props.home_upward, props.home_upward) == 0) {
        props.home_upward[1] = 1;
    }
    if (parcc_float4_dot(props.map_plane, props.map_plane) == 0) {
        props.map_plane[2] = 1;
    }
    if (props.orbit_speed[0] == 0) {
        props.orbit_speed[0] = 0.01;
    }
    if (props.orbit_speed[1] == 0) {
        props.orbit_speed[1] = 0.01;
    }
    if (props.orbit_zoom_speed == 0) {
        props.orbit_zoom_speed = 0.01;
    }

    const bool more_constrained = (int)props.map_constraint > (int)context->props.map_constraint;
    const bool orientation_changed = props.fov_orientation != context->props.fov_orientation;
    const bool viewport_resized = props.viewport_height != context->props.viewport_height ||
        props.viewport_width != context->props.viewport_width;

    context->props = props;

    if (more_constrained || orientation_changed ||
        (viewport_resized && context->props.map_constraint == PARCC_CONSTRAIN_FULL)) {
        parcc_move_with_constraints(context, context->eyepos, context->target);
    }
}

void parcc_destroy_context(parcc_context* context) { PARCC_FREE(context); }

void parcc_get_matrices(const parcc_context* context, parcc_float projection[16],
    parcc_float view[16]) {
    parcc_float gaze[3];
    parcc_float3_subtract(gaze, context->target, context->eyepos);
    parcc_float3_normalize(gaze);

    parcc_float right[3];
    parcc_float3_cross(right, gaze, context->props.home_upward);
    parcc_float3_normalize(right);

    parcc_float upward[3];
    parcc_float3_cross(upward, right, gaze);
    parcc_float3_normalize(upward);

    parcc_float16_look_at(view, context->eyepos, context->target, upward);
    const parcc_properties props = context->props;
    const parcc_float aspect = (parcc_float)props.viewport_width / props.viewport_height;
    const parcc_float fov = props.fov_degrees;
    if (context->props.fov_orientation == PARCC_HORIZONTAL) {
        parcc_float16_perspective_x(projection, fov, aspect, props.near_plane, props.far_plane);
    } else {
        parcc_float16_perspective_y(projection, fov, aspect, props.near_plane, props.far_plane);
    }
}

void parcc_get_look_at(const parcc_context* ctx, parcc_float eyepos[3], parcc_float target[3],
    parcc_float upward[3]) {
    parcc_float3_copy(eyepos, ctx->eyepos);
    parcc_float3_copy(target, ctx->target);
    if (upward) {
        parcc_float gaze[3];
        parcc_float3_subtract(gaze, ctx->target, ctx->eyepos);
        parcc_float3_normalize(gaze);

        parcc_float right[3];
        parcc_float3_cross(right, gaze, ctx->props.home_upward);
        parcc_float3_normalize(right);

        parcc_float3_cross(upward, right, gaze);
        parcc_float3_normalize(upward);
    }
}

void parcc_grab_begin(parcc_context* context, int winx, int winy, bool strafe) {
    context->grabbing = strafe ? PARCC_GRAB_STRAFE : PARCC_GRAB;

    if (context->props.mode == PARCC_MAP) {
        if (!parcc_raycast(context, winx, winy, context->grab_point_world)) {
            return;
        }
        parcc_get_ray_far(context, winx, winy, context->grab_point_far);
    }

    if (context->props.mode == PARCC_ORBIT) {
        context->grab_frame = parcc_get_current_frame(context);
        context->grab_winx = winx;
        context->grab_winy = winy;
        parcc_float3_copy(context->grab_point_pivot, context->orbit_pivot);
    }

    parcc_float3_copy(context->grab_point_eyepos, context->eyepos);
    parcc_float3_copy(context->grab_point_target, context->target);
}

void parcc_grab_update(parcc_context* context, int winx, int winy) {
    if (context->props.mode == PARCC_MAP && context->grabbing == PARCC_GRAB) {
        parcc_float u_vec[3];
        parcc_float3_subtract(u_vec, context->grab_point_world, context->grab_point_eyepos);
        const parcc_float u_len = parcc_float3_length(u_vec);

        parcc_float v_vec[3];
        parcc_float3_subtract(v_vec, context->grab_point_far, context->grab_point_world);
        const parcc_float v_len = parcc_float3_length(v_vec);

        parcc_float far_point[3];
        parcc_get_ray_far(context, winx, winy, far_point);

        parcc_float translation[3];
        parcc_float3_subtract(translation, far_point, context->grab_point_far);
        parcc_float3_scale(translation, -u_len / v_len);

        parcc_float eyepos[3];
        parcc_float3_add(eyepos, context->grab_point_eyepos, translation);

        parcc_float target[3];
        parcc_float3_add(target, context->grab_point_target, translation);

        parcc_move_with_constraints(context, eyepos, target);
    }

    if (context->props.mode == PARCC_ORBIT && context->grabbing == PARCC_GRAB) {
        parcc_frame frame = parcc_get_current_frame(context);

        const int delx = context->grab_winx - winx;
        const int dely = context->grab_winy - winy;

        const parcc_float phi = dely * context->props.orbit_speed[1];
        const parcc_float theta = delx * context->props.orbit_speed[0];

        frame.phi = context->grab_frame.phi + phi;
        frame.theta = context->grab_frame.theta + theta;
        frame.phi = PARCC_CLAMP(frame.phi, -PARCC_MAX_PHI, PARCC_MAX_PHI);

        parcc_goto_frame(context, frame);
    }

    if (context->props.mode == PARCC_ORBIT && context->grabbing == PARCC_GRAB_STRAFE) {
        parcc_float upward[3];

        parcc_float gaze[3];
        parcc_float3_subtract(gaze, context->target, context->eyepos);
        parcc_float3_normalize(gaze);

        parcc_float right[3];
        parcc_float3_cross(right, gaze, context->props.home_upward);
        parcc_float3_normalize(right);

        parcc_float3_cross(upward, right, gaze);
        parcc_float3_normalize(upward);

        const int delx = context->grab_winx - winx;
        const int dely = context->grab_winy - winy;

        const parcc_float dx = delx * 0.001;
        const parcc_float dy = dely * 0.001;

        parcc_float3_scale(right, dx);
        parcc_float3_scale(upward, dy);

        parcc_float movement[3];
        parcc_float3_add(movement, upward, right);

        parcc_float3_add(context->orbit_pivot, context->grab_point_pivot, movement);
        parcc_float3_add(context->eyepos, context->grab_point_eyepos, movement);
        parcc_float3_add(context->target, context->grab_point_target, movement);
    }
}

void parcc_zoom(parcc_context* context, int winx, int winy, parcc_float scrolldelta) {
    if (context->props.mode == PARCC_MAP) {
        parcc_float grab_point_world[3];
        if (!parcc_raycast(context, winx, winy, grab_point_world)) {
            return;
        }

        parcc_float grab_point_far[3];
        parcc_get_ray_far(context, winx, winy, grab_point_far);

        // We intentionally avoid normalizing this vector since you usually
        // want to slow down when approaching the surface.
        parcc_float u_vec[3];
        parcc_float3_subtract(u_vec, grab_point_world, context->eyepos);

        // Prevent getting stuck; this needs to be done regardless
        // of the user's min_distance setting, which is enforced in
        // parcc_move_with_constraints.
        const parcc_float zoom_speed = context->props.zoom_speed;
        if (scrolldelta > 0.0) {
            const parcc_float distance_to_surface = parcc_float3_length(u_vec);
            if (distance_to_surface < zoom_speed) {
                return;
            }
        }

        parcc_float3_scale(u_vec, scrolldelta * zoom_speed);

        parcc_float eyepos[3];
        parcc_float3_add(eyepos, context->eyepos, u_vec);

        parcc_float target[3];
        parcc_float3_add(target, context->target, u_vec);

        parcc_move_with_constraints(context, eyepos, target);
    }

    if (context->props.mode == PARCC_ORBIT) {
        parcc_float gaze[3];
        parcc_float3_subtract(gaze, context->target, context->eyepos);
        parcc_float3_normalize(gaze);
        parcc_float3_scale(gaze, context->props.orbit_zoom_speed * scrolldelta);

        parcc_float v0[3];
        parcc_float3_subtract(v0, context->orbit_pivot, context->eyepos);

        parcc_float3_add(context->eyepos, context->eyepos, gaze);
        parcc_float3_add(context->target, context->target, gaze);

        parcc_float v1[3];
        parcc_float3_subtract(v1, context->orbit_pivot, context->eyepos);

        if (parcc_float3_dot(v0, v1) < 0) {
            context->orbit_flipped = !context->orbit_flipped;
        }
    }
}

void parcc_grab_end(parcc_context* context) { context->grabbing = PARCC_GRAB_NONE; }

bool parcc_raycast(parcc_context* context, int winx, int winy, parcc_float result[3]) {
    const parcc_float width = context->props.viewport_width;
    const parcc_float height = context->props.viewport_height;
    const parcc_float fov = context->props.fov_degrees * PARCC_PI / 180.0;
    const bool vertical_fov = context->props.fov_orientation == PARCC_VERTICAL;
    const parcc_float* origin = context->eyepos;

    parcc_float gaze[3];
    parcc_float3_subtract(gaze, context->target, origin);
    parcc_float3_normalize(gaze);

    parcc_float right[3];
    parcc_float3_cross(right, gaze, context->props.home_upward);
    parcc_float3_normalize(right);

    parcc_float upward[3];
    parcc_float3_cross(upward, right, gaze);
    parcc_float3_normalize(upward);

    // Remap the grid coordinate into [-1, +1] and shift it to the pixel center.
    const parcc_float u = 2.0 * (winx + 0.5) / width - 1.0;
    const parcc_float v = 2.0 * (winy + 0.5) / height - 1.0;

    // Compute the tangent of the field-of-view angle as well as the aspect ratio.
    const parcc_float tangent = tan(fov / 2.0);
    const parcc_float aspect = width / height;

    // Adjust the gaze so it goes through the pixel of interest rather than the grid center.
    if (vertical_fov) {
        parcc_float3_scale(right, tangent * u * aspect);
        parcc_float3_scale(upward, tangent * v);
    } else {
        parcc_float3_scale(right, tangent * u);
        parcc_float3_scale(upward, tangent * v / aspect);
    }
    parcc_float3_add(gaze, gaze, right);
    parcc_float3_add(gaze, gaze, upward);
    parcc_float3_normalize(gaze);

    // Invoke the user's callback or fallback function.
    parcc_raycast_fn callback = context->props.raycast_function;
    parcc_raycast_fn fallback = parcc_raycast_plane;
    void* userdata = context->props.raycast_userdata;
    if (!callback) {
        callback = fallback;
        userdata = context;
    }

    // If the ray misses, then try the fallback function.
    parcc_float t;
    if (!callback(origin, gaze, &t, userdata)) {
        if (callback == fallback) {
            return false;
        }
        if (!fallback(origin, gaze, &t, context)) {
            return false;
        }
    }

    parcc_float3_scale(gaze, t);
    parcc_float3_add(result, origin, gaze);
    return true;
}

parcc_frame parcc_get_current_frame(const parcc_context* context) {
    parcc_frame frame;
    frame.mode = context->props.mode;

    if (context->props.mode == PARCC_MAP) {
        const parcc_float* origin = context->eyepos;
        const parcc_float* upward = context->props.home_upward;

        parcc_float direction[3];
        parcc_float3_subtract(direction, context->target, origin);
        parcc_float3_normalize(direction);

        parcc_float distance;
        parcc_raycast_plane(origin, direction, &distance, (void*)context);

        const parcc_float fov = context->props.fov_degrees * PARCC_PI / 180.0;
        const parcc_float half_extent = distance * tan(fov / 2);

        parcc_float target[3];
        parcc_float3_scale(direction, distance);
        parcc_float3_add(target, origin, direction);

        // Compute the tangent frame defined by the map_plane normal and the home_upward vector.
        parcc_float uvec[3];
        parcc_float vvec[3];
        parcc_float target_to_eye[3];

        parcc_float3_copy(target_to_eye, context->props.map_plane);
        parcc_float3_cross(uvec, upward, target_to_eye);
        parcc_float3_cross(vvec, target_to_eye, uvec);

        parcc_float3_subtract(target, target, context->props.home_target);

        frame.extent = half_extent * 2;
        frame.center[0] = parcc_float3_dot(uvec, target);
        frame.center[1] = parcc_float3_dot(vvec, target);
    }

    if (context->props.mode == PARCC_ORBIT) {
        parcc_float pivot_to_eye[3];
        parcc_float3_subtract(pivot_to_eye, context->eyepos, context->orbit_pivot);
        const parcc_float d = parcc_float3_length(pivot_to_eye);
        const parcc_float x = pivot_to_eye[0] / d;
        const parcc_float y = pivot_to_eye[1] / d;
        const parcc_float z = pivot_to_eye[2] / d;
        frame.phi = asin(y);
        frame.theta = atan2(x, z);
        frame.pivot_distance = context->orbit_flipped ? -d : d;
        parcc_float3_copy(frame.pivot, context->orbit_pivot);
    }

    return frame;
}

parcc_frame parcc_get_home_frame(const parcc_context* context) {
    const parcc_float width = context->props.viewport_width;
    const parcc_float height = context->props.viewport_height;
    const parcc_float aspect = width / height;

    parcc_frame frame;
    frame.mode = context->props.mode;

    if (frame.mode == PARCC_MAP) {
        const parcc_float map_width = context->props.map_extent[0] / 2;
        const parcc_float map_height = context->props.map_extent[1] / 2;
        const bool horiz = context->props.fov_orientation == PARCC_HORIZONTAL;
        frame.extent = horiz ? context->props.map_extent[0] : context->props.map_extent[1];
        frame.center[0] = 0;
        frame.center[1] = 0;
        if (context->props.map_constraint != PARCC_CONSTRAIN_FULL) {
            return frame;
        }
        if (horiz) {
            parcc_float vp_width = frame.extent / 2;
            parcc_float vp_height = vp_width / aspect;
            if (map_height < vp_height) {
                frame.extent = 2 * map_height * aspect;
            }
        } else {
            parcc_float vp_height = frame.extent / 2;
            parcc_float vp_width = vp_height * aspect;
            if (map_width < vp_width) {
                frame.extent = 2 * map_width / aspect;
            }
        }
    }

    if (frame.mode == PARCC_ORBIT) {
        frame.theta = frame.phi = 0;
        parcc_float3_copy(frame.pivot, context->props.home_target);
        frame.pivot_distance = parcc_float3_length(context->props.home_vector);
    }

    return frame;
}

void parcc_goto_frame(parcc_context* context, parcc_frame frame) {
    if (context->props.mode == PARCC_MAP) {
        const parcc_float* upward = context->props.home_upward;
        const parcc_float half_extent = frame.extent / 2.0;
        const parcc_float fov = context->props.fov_degrees * PARCC_PI / 180.0;
        const parcc_float distance = half_extent / tan(fov / 2);

        // Compute the tangent frame defined by the map_plane normal and the home_upward vector.
        parcc_float uvec[3];
        parcc_float vvec[3];
        parcc_float target_to_eye[3];

        parcc_float3_copy(target_to_eye, context->props.map_plane);
        parcc_float3_cross(uvec, upward, target_to_eye);
        parcc_float3_cross(vvec, target_to_eye, uvec);

        // Scale the U and V components by the frame coordinate.
        parcc_float3_scale(uvec, frame.center[0]);
        parcc_float3_scale(vvec, frame.center[1]);

        // Obtain the new target position by adding U and V to home_target.
        parcc_float3_copy(context->target, context->props.home_target);
        parcc_float3_add(context->target, context->target, uvec);
        parcc_float3_add(context->target, context->target, vvec);

        // Obtain the new eye position by adding the scaled plane normal to the new target
        // position.
        parcc_float3_scale(target_to_eye, distance);
        parcc_float3_add(context->eyepos, context->target, target_to_eye);
    }

    if (context->props.mode == PARCC_ORBIT) {
        parcc_float3_copy(context->orbit_pivot, frame.pivot);
        const parcc_float x = sin(frame.theta) * cos(frame.phi);
        const parcc_float y = sin(frame.phi);
        const parcc_float z = cos(frame.theta) * cos(frame.phi);
        parcc_float3_set(context->eyepos, x, y, z);
        parcc_float3_scale(context->eyepos, fabs(frame.pivot_distance));
        parcc_float3_add(context->eyepos, context->eyepos, context->orbit_pivot);

        context->orbit_flipped = frame.pivot_distance < 0;

        parcc_float3_set(context->target, x, y, z);
        parcc_float3_scale(context->target, context->orbit_flipped ? 1.0 : -1.0);
        parcc_float3_add(context->target, context->target, context->eyepos);
    }
}

parcc_frame parcc_interpolate_frames(parcc_frame a, parcc_frame b, double t) {
    const double rho = sqrt(2.0);
    const double rho2 = 2, rho4 = 4;
    const double ux0 = a.center[0], uy0 = a.center[1], w0 = a.extent;
    const double ux1 = b.center[0], uy1 = b.center[1], w1 = b.extent;
    const double dx = ux1 - ux0, dy = uy1 - uy0, d2 = dx * dx + dy * dy, d1 = sqrt(d2);
    const double b0 = (w1 * w1 - w0 * w0 + rho4 * d2) / (2.0 * w0 * rho2 * d1);
    const double b1 = (w1 * w1 - w0 * w0 - rho4 * d2) / (2.0 * w1 * rho2 * d1);
    const double r0 = log(sqrt(b0 * b0 + 1.0) - b0);
    const double r1 = log(sqrt(b1 * b1 + 1) - b1);
    const double dr = r1 - r0;
    const int valid_dr = (dr == dr) && dr != 0;
    const double S = (valid_dr ? dr : log(w1 / w0)) / rho;
    parcc_frame frame;
    const double s = t * S;
    if (valid_dr) {
        const double coshr0 = cosh(r0);
        const double u = w0 / (rho2 * d1) * (coshr0 * tanh(rho * s + r0) - sinh(r0));
        frame.center[0] = ux0 + u * dx;
        frame.center[1] = uy0 + u * dy;
        frame.extent = w0 * coshr0 / cosh(rho * s + r0);
        return frame;
    }
    frame.center[0] = ux0 + t * dx;
    frame.center[1] = uy0 + t * dy;
    frame.extent = w0 * exp(rho * s);
    return frame;
}

double parcc_get_interpolation_duration(parcc_frame a, parcc_frame b) {
    const double rho = sqrt(2.0);
    const double rho2 = 2, rho4 = 4;
    const double ux0 = a.center[0], uy0 = a.center[1], w0 = a.extent;
    const double ux1 = b.center[0], uy1 = b.center[1], w1 = b.extent;
    const double dx = ux1 - ux0, dy = uy1 - uy0, d2 = dx * dx + dy * dy, d1 = sqrt(d2);
    const double b0 = (w1 * w1 - w0 * w0 + rho4 * d2) / (2.0 * w0 * rho2 * d1);
    const double b1 = (w1 * w1 - w0 * w0 - rho4 * d2) / (2.0 * w1 * rho2 * d1);
    const double r0 = log(sqrt(b0 * b0 + 1.0) - b0);
    const double r1 = log(sqrt(b1 * b1 + 1) - b1);
    const double dr = r1 - r0;
    const int valid_dr = (dr == dr) && dr != 0;
    const double S = (valid_dr ? dr : log(w1 / w0)) / rho;
    return fabs(S);
}

static bool parcc_raycast_plane(const parcc_float origin[3], const parcc_float dir[3],
    parcc_float* t, void* userdata) {
    parcc_context* context = (parcc_context*)userdata;
    const parcc_float* plane = context->props.map_plane;
    parcc_float n[3] = {plane[0], plane[1], plane[2]};
    parcc_float p0[3] = {plane[0], plane[1], plane[2]};
    parcc_float3_scale(p0, plane[3]);
    const parcc_float denom = -parcc_float3_dot(n, dir);
    if (denom > 1e-6) {
        parcc_float p0l0[3];
        parcc_float3_subtract(p0l0, p0, origin);
        *t = parcc_float3_dot(p0l0, n) / -denom;
        return *t >= 0;
    }
    return false;
}

// Finds the point on the frustum's far plane that a pick ray intersects.
static void parcc_get_ray_far(parcc_context* context, int winx, int winy, parcc_float result[3]) {
    const parcc_float width = context->props.viewport_width;
    const parcc_float height = context->props.viewport_height;
    const parcc_float fov = context->props.fov_degrees * PARCC_PI / 180.0;
    const bool vertical_fov = context->props.fov_orientation == PARCC_VERTICAL;
    const parcc_float* origin = context->eyepos;

    parcc_float gaze[3];
    parcc_float3_subtract(gaze, context->target, origin);
    parcc_float3_normalize(gaze);

    parcc_float right[3];
    parcc_float3_cross(right, gaze, context->props.home_upward);
    parcc_float3_normalize(right);

    parcc_float upward[3];
    parcc_float3_cross(upward, right, gaze);
    parcc_float3_normalize(upward);

    // Remap the grid coordinate into [-1, +1] and shift it to the pixel center.
    const parcc_float u = 2.0 * (winx + 0.5) / width - 1.0;
    const parcc_float v = 2.0 * (winy + 0.5) / height - 1.0;

    // Compute the tangent of the field-of-view angle as well as the aspect ratio.
    const parcc_float tangent = tan(fov / 2.0);
    const parcc_float aspect = width / height;

    // Adjust the gaze so it goes through the pixel of interest rather than the grid center.
    if (vertical_fov) {
        parcc_float3_scale(right, tangent * u * aspect);
        parcc_float3_scale(upward, tangent * v);
    } else {
        parcc_float3_scale(right, tangent * u);
        parcc_float3_scale(upward, tangent * v / aspect);
    }
    parcc_float3_add(gaze, gaze, right);
    parcc_float3_add(gaze, gaze, upward);
    parcc_float3_scale(gaze, context->props.far_plane);
    parcc_float3_add(result, origin, gaze);
}

static void parcc_move_with_constraints(parcc_context* context, const parcc_float eyepos[3],
    const parcc_float target[3]) {
    const parcc_constraint constraint = context->props.map_constraint;
    const parcc_float width = context->props.viewport_width;
    const parcc_float height = context->props.viewport_height;
    const parcc_float aspect = width / height;
    const parcc_float map_width = context->props.map_extent[0] / 2;
    const parcc_float map_height = context->props.map_extent[1] / 2;
    const parcc_frame home = parcc_get_home_frame(context);
    const parcc_frame previous_frame = parcc_get_current_frame(context);
    const parcc_float fov = context->props.fov_degrees * PARCC_PI / 180.0;
    const parcc_float min_extent = 2.0 * context->props.map_min_distance * tan(fov / 2);

    parcc_float3_copy(context->eyepos, eyepos);
    parcc_float3_copy(context->target, target);

    parcc_frame frame = parcc_get_current_frame(context);

    if (frame.extent < min_extent) {
        frame.extent = min_extent;
        frame.center[0] = previous_frame.center[0];
        frame.center[1] = previous_frame.center[1];
    }

    if (constraint == PARCC_CONSTRAIN_NONE) {
        parcc_goto_frame(context, frame);
        return;
    }

    parcc_float x = frame.center[0];
    parcc_float y = frame.center[1];

    if (context->props.fov_orientation == PARCC_HORIZONTAL) {
        parcc_float vp_width = frame.extent / 2;
        parcc_float vp_height = vp_width / aspect;
        if (map_width < vp_width) {
            frame.extent = home.extent;
            vp_width = frame.extent / 2;
            vp_height = vp_width / aspect;
            x = 0;
            y = previous_frame.center[1];
        }
        x = PARCC_CLAMP(x, -map_width + vp_width, map_width - vp_width);
        if (map_height < vp_height) {
            if (context->props.map_constraint == PARCC_CONSTRAIN_FULL) {
                frame.extent = 2 * map_height * aspect;
                vp_width = frame.extent / 2;
                vp_height = vp_width / aspect;
                x = previous_frame.center[0];
                x = PARCC_CLAMP(x, -map_width + vp_width, map_width - vp_width);
                y = PARCC_CLAMP(y, -map_height + vp_height, map_height - vp_height);
            } else {
                y = PARCC_CLAMP(y, -vp_height + map_height, vp_height - map_height);
            }
        } else {
            y = PARCC_CLAMP(y, -map_height + vp_height, map_height - vp_height);
        }
    } else {
        parcc_float vp_height = frame.extent / 2;
        parcc_float vp_width = vp_height * aspect;
        if (map_height < vp_height) {
            frame.extent = home.extent;
            vp_height = frame.extent / 2;
            vp_width = vp_height * aspect;
            y = 0;
            x = previous_frame.center[0];
        }
        y = PARCC_CLAMP(y, -map_height + vp_height, map_height - vp_height);
        if (map_width < vp_width) {
            if (context->props.map_constraint == PARCC_CONSTRAIN_FULL) {
                frame.extent = 2 * map_width / aspect;
                vp_height = frame.extent / 2;
                vp_width = vp_height * aspect;
                y = previous_frame.center[1];
                y = PARCC_CLAMP(y, -map_height + vp_height, map_height - vp_height);
                x = PARCC_CLAMP(x, -map_width + vp_width, map_width - vp_width);
            } else {
                x = PARCC_CLAMP(x, -vp_width + map_width, vp_width - map_width);
            }
        } else {
            x = PARCC_CLAMP(x, -map_width + vp_width, map_width - vp_width);
        }
    }

    frame.center[0] = x;
    frame.center[1] = y;
    parcc_goto_frame(context, frame);
}

#endif  // PAR_CAMERA_CONTROL_IMPLEMENTATION
#endif  // PAR_CAMERA_CONTROL_H

// par_camera_control is distributed under the MIT license:
//
// Copyright (c) 2019 Philip Rideout
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
