/**
 * Path tracer data structures, but set up for sending to an
 * OpenCL device (padding and all).
 *
 * (c) L. Diener 2013
 */

 #ifndef __DATA_STRUCTURES_H__
 #define __DATA_STRUCTURES_H__
 
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
	return (vec3_t){x, y, z, 0.0f}; 
}

// Data for a round sphere
// Padded to align on 4 floats.
typedef struct sphere_object_data {
    vec3_t position;
    float radius;
	float pad_a;
	float pad_b;
	float pad_c;
} sphere_object_data;

// Data for a brdf that really doesn't make much of an effort
// Padded to align on 4 floats.
typedef struct simple_brdf_data {
    vec3_t emission;
    vec3_t albedo;
    float reflectivity;
    float specularity;
	float pad_a;
	float pad_b;
} simple_brdf_data;

// A thing that can be displayed
typedef struct object {
    sphere_object_data object_data; // Data for object
    simple_brdf_data brdf_data; // Data for brdf
} object;

#endif