/**
 * Simple OpenCL path tracing kernel.
 *
 * (c) L. Diener, 2013
 */

// Objects in scene
#define OBJECT_COUNT 13

// Maximum hit stack depth - needs to be big to avoid bias!
// How much is enough depends on termination probability. 
// For 0.95, 1024 is plenty.
#define MAX_DEPTH 1024

// Likelihood of continuing a ray. Big numbers means slower
// sampling but usually more variance reduction per sample.
// The ideal value depends on your scene. 0.95 is a good
// choice. Make sure the hit stack is big enough!
#define CONTINUATION_PROBABILITY 0.95f;

// Image size. Even numbers, please.
#define IMAGE_WIDTH 400
#define IMAGE_HEIGHT 400

// Samples per pixel
#define SAMPLES 30000

// LFSR data
typedef struct {
	ulong a;
	ulong b;
	ulong c;
} random_state;

// A simple ray
typedef struct ray {
    float3 origin; // From where
    float3 direction; // To where
    float3 radiance; // And what
    float terminate; // Until when
} ray;

// A hit, or maybe not
typedef struct hit {
    int object; // Hit what (via id in scene)
    float3 position; // Hit where
    float3 normal; // Points away
    float3 s; // Tangent space
    float3 t; // Also tangent space
    ray hitter; // Ray that went in
    int hit; // 1: yep, 0: nope
    float distance; // How far did we have to walk the ray
} hit;

// Data for a brdf that really doesn't make much of an effort
typedef struct simple_brdf_data {
    float3 emission;
    float3 albedo;
    float reflectivity;
    float specularity;
} simple_brdf_data;

// Data for a round sphere
typedef struct sphere_object_data {
    float3 position;
    float radius;
} sphere_object_data;

// A thing that can be displayed
typedef struct object {
    sphere_object_data object_data; // Data for object
    simple_brdf_data brdf_data; // Data for brdf
} object;

// Sampling, uniform-unit
inline float sample_unit(random_state *r) {
	ulong old = r->b;
	r->b = r->a * 1103515245 + 12345;
	r->a = (~old ^ (r->b >> 3)) - r->c++;
    return (float)(r->b & 4294967295) / 4294967295.0f;
}

// LFSR seeder
void seed_random(random_state *r, uint seed) {
	r->a = seed;
	r->b = 0;
	r->c = 362436;
	for( int i = 0; i < 1000; i++ ) {
		sample_unit( r );
	}
}

// Sampling, uniform from a hemisphere oriented towards a vector
inline float3 sample_hemisphere_uniform(random_state *r, float3 towards) {
    float u = sample_unit(r) * 2.0f - 1.0f;
    float v = sample_unit(r) * 2.0f * 3.14159265359f;
	float w = sqrt(1.0f - pow(u, 2.0f));
    float3 sample = (float3)(w * cos(v), w * sin(v), u);
    return dot(sample, towards) < 0 ? -sample : sample;
}

// Sampling, from a hemisphere oriented towards a vector, weighted to center
float3 sample_hemisphere_weighted(random_state *r, float3 towards, float weight) {
    float3 sample = sample_hemisphere_uniform(r, towards);
    return normalize(sample + towards * weight);
}

// Shove point in direction, to avoid speckles
inline float3 displace(float3 in, float3 towards) {
    return in + towards * 0.002f;
}

// It's a sphere. It's round
inline hit sphere_intersect(__constant sphere_object_data* object_data, ray ray_in) {    
    hit hit_in;

    // Begin intersection test
    float3 x = ray_in.origin - object_data->position;
    float a = dot(ray_in.direction, ray_in.direction);
    float b = 2.0f * dot( ray_in.direction, x );
    float c = dot(x, x) - (object_data->radius * object_data->radius);
    float discriminant = b * b - 4.0f * a * c;

    // Condition 1: Discriminant negative
    if( discriminant < 0 ) {
        hit_in.hit = 0;
        return hit_in;
    }
    
    float discsqrt = sqrt(discriminant);
    float q;
    if (b < 0) {
        q = (-b - discsqrt) / 2.0f;
    }
    else {
        q = (-b + discsqrt) / 2.0f;
    }

    // Make t1 the bigger distance
    float t0 = q / a;
    float t1 = c / q;
    if(t0 > t1) {
        float temp = t0;
        t0 = t1;
        t1 = temp;
    }

    // Condition 2: Bigger distance is still negative
    if (t1 < 0) {
        hit_in.hit = 0;
        return hit_in;
    }

    // Seems like we hit, figure out which hit is closer
    float intersection_distance;
    if (t0 < 0) {
        intersection_distance = t1;
    }
    else {
        intersection_distance = t0;
    }

    // Fill up the hit struct
    hit_in.position = ray_in.origin + ray_in.direction * intersection_distance;
    hit_in.normal = normalize(hit_in.position - object_data->position);
    if(hit_in.normal.x <= 0.1f && hit_in.normal.y <= 0.1f) {
        hit_in.s = normalize((float3)(hit_in.normal.x, hit_in.normal.y, 0));
    }
    else {
        hit_in.s = normalize((float3)(hit_in.normal.x, 0, hit_in.normal.z));
    }
    hit_in.t = normalize(cross(hit_in.normal, hit_in.s));
    hit_in.hitter = ray_in;
    hit_in.hit = 1;
    hit_in.distance = intersection_distance;
    
    return hit_in;
}

// Phenomenological and pulled out of ass, mostly. Lambertian, I hope.
inline ray simple_brdf_reflect(random_state *r, __constant simple_brdf_data* brdf_data, hit hit_in) {
    ray ray_out;
    
    if(sample_unit(r) > brdf_data->reflectivity) {
        // Next ray is diffuse
        float3 out_dir = sample_hemisphere_uniform(r, hit_in.normal);
        ray_out.direction = normalize(out_dir);
        ray_out.origin = displace(hit_in.position, ray_out.direction);
    }
    else {
        // Next ray is specular
        float3 axis = sample_hemisphere_weighted(r, hit_in.normal, brdf_data->specularity);
        ray_out.direction = normalize(hit_in.hitter.direction - axis * dot(hit_in.hitter.direction, axis) * 2.0f);
        ray_out.origin = displace(hit_in.position, ray_out.direction);
    }
    
    ray_out.terminate = hit_in.hitter.terminate * CONTINUATION_PROBABILITY;
    return ray_out;
}

// The Cool Thing is that the lambertian term is all in the sampling so 
// the actual shading is very simple.
inline float3 simple_brdf_shade(__constant simple_brdf_data* brdf_data, float3 radiance_out) {
    float3 radiance = radiance_out;
    radiance = radiance + brdf_data->emission;
    radiance = radiance * brdf_data->albedo;
    return radiance;
}

// Tracing. Put a ray in, get radiance along ray out. Sample a lot for better results.
inline float3 shade(random_state *r, ray ray_in, __constant object* scene, float* ray_probabilities, int* hit_objects) {
    // Iterate over rays until russian roulette decides it's over or we run out of stack
    int depth = 0;
    ray current_ray = ray_in;
	hit best_hit;
	ray_probabilities[0] = ray_in.terminate;
    while(depth < MAX_DEPTH && sample_unit(r) <= current_ray.terminate) {   
        best_hit.hit = 0;
        best_hit.distance = FLT_MAX;
        best_hit.object = 0;
        
        for(int i = 0; i < OBJECT_COUNT; i++) {
            hit current_hit = sphere_intersect(&scene[i].object_data, current_ray);
            if(current_hit.hit == 1 && current_hit.distance < best_hit.distance) {
                best_hit = current_hit;
                best_hit.object = i;
            }
        }
        hit_objects[depth] = best_hit.object;
		depth++;
		ray_probabilities[depth] = current_ray.terminate;		
        current_ray = simple_brdf_reflect(r, &scene[best_hit.object].brdf_data, best_hit);
    }
    
    // Go backwards from the end ray and shade
	float3 radiance = (float3)(0, 0, 0);
    for(int backwards_ray = depth-1; backwards_ray >= 0; backwards_ray--) {
        radiance = simple_brdf_shade(
            &scene[hit_objects[backwards_ray]].brdf_data,
            radiance
        );
		radiance = radiance * ray_probabilities[backwards_ray];
    }

    return radiance;
}

__kernel void path_tracer(__global float* r, __global float* g, __global float* b, __constant object* scene) {
	// Coordinates
	unsigned int i = get_global_id(0);
	int coord_x = (i % IMAGE_WIDTH) - IMAGE_WIDTH / 2;
	int coord_y = -(floor((float)i / (float)IMAGE_WIDTH) - IMAGE_HEIGHT / 2);
	
	// Randomness
	random_state randstate;
	seed_random(&randstate, i);
	
	// Prepare hit stack
	float ray_probabilities[MAX_DEPTH];
    int hit_objects[MAX_DEPTH];
	
	// Trace
	float3 result = (float3)(0, 0, 0);
	ray trace_ray;
	for(int i = 0; i < SAMPLES; i++) {
		// Generate ray with fuzzed pixel coordinates for AA.
		trace_ray.origin = (float3)(0, 0, 0);
		trace_ray.direction = normalize((float3)(
			(coord_x + sample_unit(&randstate) - 0.5f) * (800.0f / (float)IMAGE_WIDTH),
			(coord_y + sample_unit(&randstate) - 0.5f) * (800.0f / (float)IMAGE_HEIGHT),
			-800.0f
		));
		trace_ray.terminate = 1.0f;	
		result += shade(&randstate, trace_ray, scene, ray_probabilities, hit_objects);
	}
	
	r[i] = result.x;
	g[i] = result.y;
	b[i] = result.z;
}
