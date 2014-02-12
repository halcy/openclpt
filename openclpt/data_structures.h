/**
 * Path tracer data structures, but set up for sending to an
 * OpenCL device (padding and all).
 *
 * (c) L. Diener 2013
 */

 #ifndef __DATA_STRUCTURES_H__
 #define __DATA_STRUCTURES_H__
 
#include <math.h>

// A three-float vector, represented as a four-float
// vector because OpenCL internals.
typedef struct vec3_t {
	float x;
	float y;
	float z;
	float w;
} vec3_t;
 
 // A "constructor"
static inline vec3_t vec3(float x, float y, float z) {
	vec3_t v;
	v.x = x;
	v.y = y;
	v.z = z;
	v.w = 0.0f;
	return v; 
}

// Some math functions
static inline vec3_t add(vec3_t a, vec3_t b) {
	return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline vec3_t sub(vec3_t a, vec3_t b) {
	return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline vec3_t mul(vec3_t a, float b) {
	return vec3(a.x * b, a.y * b, a.z * b);
}

static inline vec3_t div(vec3_t a, float b) {
	return vec3(a.x / b, a.y / b, a.z / b);
}

static inline float dot(vec3_t a, vec3_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline float length(vec3_t v) {
	return sqrt(dot(v, v));
}

static inline vec3_t norm(vec3_t v) {
	return div(v, length(v));
}

static inline vec3_t cross(vec3_t a, vec3_t b) {
	return vec3(
		a.y * b.z - a.z * b.y, 
		a.z * b.x - a.x * b.z, 
		a.x * b.y - a.y * b.x
	);
}

static inline vec3_t project(vec3_t a, vec3_t n) {
	return mul(n, dot(a, n));
}

// Data for an AABB
// Vertices should be ordered smallest -> biggest
// Use the "aabb" function to make them to ensure this
typedef struct aabb_data {
	vec3_t vertices[8];
} aabb_data;

static inline aabb_data aabb(vec3_t center, vec3_t radius) {
	aabb_data ret;
	ret.vertices[0] = add(center, vec3(-radius.x, -radius.y, -radius.z));
	ret.vertices[1] = add(center, vec3(-radius.x, -radius.y,  radius.z));
	ret.vertices[2] = add(center, vec3(-radius.x,  radius.y, -radius.z));
	ret.vertices[3] = add(center, vec3(-radius.x,  radius.y,  radius.z));
	ret.vertices[4] = add(center, vec3( radius.x, -radius.y, -radius.z));
	ret.vertices[5] = add(center, vec3( radius.x, -radius.y,  radius.z));
	ret.vertices[6] = add(center, vec3( radius.x,  radius.y, -radius.z));
	ret.vertices[7] = add(center, vec3( radius.x,  radius.y,  radius.z));
	return ret;
}


// Data for a triangle
// Padded to align on 4 floats.
typedef struct triangle_object_data {
	vec3_t p1;
	vec3_t p2;
	vec3_t p3;
	vec3_t n1;
	vec3_t n2;
	vec3_t n3;
} triangle_object_data;

// Data for a brdf that really doesn't make much of an effort
// Padded to align on 4 floats.
typedef struct simple_brdf_data {
    vec3_t emission;
    vec3_t albedo;
    float reflectivity;
	float transparency;	
    float specularity;
	float pad_b;
} simple_brdf_data;

// A thing that can be displayed
// SUPER DENORMALIZED for easy access
typedef struct object {
    triangle_object_data object_data; // Data for object
    simple_brdf_data brdf_data; // Data for brdf
} object;

#endif