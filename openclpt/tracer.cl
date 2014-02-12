/**
 * Simple OpenCL path tracing kernel.
 *
 * (c) L. Diener, 2014
 */

// Maximum hit stack depth - needs to be big to avoid bias!
// How much is enough depends on termination probability. 
// For 0.95, 1024 is plenty.
#define MAX_DEPTH 1024

// A small number
#define EPSILON 0.0001f 

// Likelihood of continuing a ray. Big numbers means slower
// sampling but usually more variance reduction per sample.
// The ideal value depends on your scene. 0.95 is a good
// choice. Make sure the hit stack is big enough!
#define CONTINUATION_PROBABILITY 0.95f

// Samples per pixel per call
#define NORMDEPTH_SAMPLES 50

// Grid stuff
#define GRID_RES_HALF (0.5f * GRID_RES)
#define GRID_RADIUS_2 (2 * GRID_RADIUS)
#define GRID_MAX (GRID_RADIUS * GRID_RES - GRID_RES_HALF)
#define GRID_ITEM(x) ((int)((x) / GRID_RES + GRID_RADIUS))
#define GRID_CELL(x, y, z) ((int)((x) + (y) * GRID_RADIUS_2 + (z) * GRID_RADIUS_2 * GRID_RADIUS_2))
#define GRID_ID(pos) (GRID_CELL(GRID_ITEM(pos.x), GRID_ITEM(pos.y), GRID_ITEM(pos.z)))
#define GRID_X(x) ((float)((x) % GRID_RADIUS_2) * GRID_RES - GRID_MAX)
#define GRID_Y(x) ((float)((x / GRID_RADIUS_2) % GRID_RADIUS_2) * GRID_RES - GRID_MAX)
#define GRID_Z(x) ((float)((x / (GRID_RADIUS_2 * GRID_RADIUS_2))) * GRID_RES - GRID_MAX)
#define GRID_POS(x) ((float3)(GRID_X(x), GRID_Y(x), GRID_Z(x)))

// Sampler type for reading image data
__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;

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
    ray hitter; // Ray that went in
    int hit; // 1: yep, 0: nope
    float distance; // How far did we have to walk the ray
} hit;

// Data for a brdf that really doesn't make much of an effort
typedef struct simple_brdf_data {
    float3 emission;
    float3 albedo;
    float reflectivity;
    float transparency;	
    float specularity;
} simple_brdf_data;

// Data for a triangle
// Padded to align on 4 floats.
typedef struct triangle_object_data {
	float3 p1;
	float3 p2;
	float3 p3;
	float3 n1;
	float3 n2;
	float3 n3;
} triangle_object_data;

// A thing that can be displayed
// Super denormalized for easy access
typedef struct object {
    triangle_object_data object_data; // Data for object
    simple_brdf_data brdf_data; // Data for brdf
} object;

// Sampling, uniform-unit [0, 1]
inline float sample_unit(random_state *r) {
	ulong old = r->b;
	r->b = r->a * 1103515245 + 12345;
	r->a = (~old ^ (r->b >> 3)) - r->c++;
    return (float)(r->b & 4294967295) / 4294967295.0f;
}

// LFSR seeder
void seed_random(random_state *r, const uint seed) {
	r->a = seed;
	r->b = 0;
	r->c = 362436;
	for( int i = 0; i < 1000; i++ ) {
		sample_unit( r );
	}
}

// Sampling, uniform from a hemisphere oriented towards a vector
inline float3 sample_hemisphere_uniform(random_state *r, const float3 towards) {
    float u = sample_unit(r) * 2.0f - 1.0f;
    float v = sample_unit(r) * 2.0f * 3.14159265359f;
	float w = sqrt(1.0f - pow(u, 2.0f));
    float3 sample = (float3)(w * cos(v), w * sin(v), u);
    return dot(sample, towards) < 0 ? -sample : sample;
}

// Sampling, from a hemisphere oriented towards a vector, weighted to center
float3 sample_hemisphere_weighted(random_state *r, float3 towards, const float weight) {
    float3 sample = sample_hemisphere_uniform(r, towards);
    return normalize(sample + towards * weight);
}

// Shove point in direction, to avoid speckles
inline float3 displace(const float3 in, const float3 towards) {
    return in + towards * EPSILON * 2.0f;
}

// Triangle intersection with barycentric normal interpolation
inline hit triangle_intersect(const __global triangle_object_data* object_data, const ray ray_in) {    
    hit hit_in;
	
	// Vertices
	float3 A = object_data->p1;
	float3 B = object_data->p2;
	float3 C = object_data->p3;

	// Figure out triangle plane
	float3 triangle_ab = B - A;
	float3 triangle_ac = C - A;
	float3 triangle_nn = cross(triangle_ab, triangle_ac);
	float3 triangle_n = normalize(triangle_nn);
	float triangle_support = dot(A, triangle_n);

	// Compute intersection distance, bail if infinite or negative
	float intersection_det = dot(triangle_n, ray_in.direction);
	if(fabs(intersection_det) <= EPSILON) {
		hit_in.hit = 0;
        return hit_in;
	}
	float intersection_dist = (triangle_support - dot(triangle_n, ray_in.origin)) / intersection_det;
	if(intersection_dist <= 0.0f) {
		hit_in.hit = 0;
        return hit_in;
	}

	// Compute intersection point
	float3 Q = ray_in.origin + ray_in.direction * intersection_dist;

	// Test inside-ness
	float3 triangle_bc = C - B;
	float3 triangle_ca = A - C;
	float3 triangle_aq = Q - A;
	float3 triangle_bq = Q - B;
	float3 triangle_cq = Q - C;

	float baryA = dot(cross(triangle_bc, triangle_bq), triangle_n);
	float baryB = dot(cross(triangle_ca, triangle_cq), triangle_n);
	float baryC = dot(cross(triangle_ab, triangle_aq), triangle_n);

	if(baryA < 0.0f || baryB < 0.0f || baryC < 0.0f) {
		hit_in.hit = 0;
        return hit_in;
	}

	// Perform barycentric interpolation of normals
	float triangle_den = dot(triangle_nn, triangle_n);
	baryA /= triangle_den;
	baryB /= triangle_den;
	baryC /= triangle_den;

	float3 N = normalize(
		object_data->n1 * baryA + 
		object_data->n2 * baryB + 
		object_data->n3 * baryC
	);

    // Fill up the hit struct
    hit_in.position = Q;
    hit_in.normal = N;
    hit_in.hitter = ray_in;
    hit_in.hit = 1;
    hit_in.distance = intersection_dist;
    
    return hit_in;
}
// A helper function: Given normal direction d, normal n and indices r1 (inside) and r2 (outside),
// return a refraction direction
float3 simple_brdf_refract(random_state* rstate, const float3 d, const float3 n, const float r1, const float r2) {
	float rr1;
	float rr2;
	float3 rn;
	
	// Figure out whether we're going in or out and flip things accordingly
	float theta1 = dot(d, n);
	if(theta1 > 0 ) {
		rr1 = r1;
		rr2 = r2;
		rn = -n;
	}
	else {
		rr1 = r2;
		rr2 = r1;
		rn = n;
		theta1 = -theta1;
	}
	
	// Figure out whether we have total internal reflection
	float r = rr1 / rr2;
	float theta2 = sqrt(1.0f - r * r * (1.0f - theta1 * theta1));
	if(theta2 < 0) {
		return normalize(d - rn * dot(d, rn) * 2.0f);
	}
	  
	// Figure out what the Fresnel equations say about what happens next
	float rs = (rr1 * theta1 - rr2 * theta2) / (rr1 * theta1 + rr2 * theta2);
	rs = rs * rs;
	float rp = (rr1 * theta2 - rr2 * theta1) / (rr1 * theta2 + rr2 * theta1);
	rp = rp * rp;
	float rr = (rs + rp) / 2.0f;
	
	// Choose to either refract or reflect, based on fresnel coefficient
	if(sample_unit(rstate) > rr) {
		// Refract
		return normalize(d * r - rn * (r * theta1 + sqrt(theta2)));
	}
	else {
		// Reflect
		return normalize(d - rn * dot(d, rn) * 2.0f);
	}
}

// Phenomenological and pulled out of ass, mostly. Lambertian, I hope.
inline ray simple_brdf_reflect(random_state* r, const __global simple_brdf_data* brdf_data, const hit hit_in) {
    ray ray_out;
    
    if(sample_unit(r) > brdf_data->reflectivity) {
        // Next ray is diffuse
        ray_out.direction = normalize(sample_hemisphere_uniform(r, hit_in.normal));
    }
    else {
        // Next ray is specular or transmission TODO FIXME reject samples that end up going "in" (unless that's fine)
		// Also see to it that handling of refractive indices is correct
        float3 axis = sample_hemisphere_weighted(r, hit_in.normal, brdf_data->specularity);
		if(sample_unit(r) > brdf_data->transparency) {
			// Reflect
			ray_out.direction = normalize(hit_in.hitter.direction - axis * dot(hit_in.hitter.direction, axis) * 2.0f);
		} 
		else {
			// Refract
			ray_out.direction = simple_brdf_refract(r, hit_in.hitter.direction, axis, 1.5f, 1.0f);
		}
    }
    
	ray_out.origin = displace(hit_in.position, ray_out.direction);
    ray_out.terminate = hit_in.hitter.terminate * CONTINUATION_PROBABILITY;
    return ray_out;
}

// The Cool Thing is that the lambertian term and whatnot is all in the sampling so 
// the actual shading is very simple.
inline float3 simple_brdf_shade(const __global simple_brdf_data* brdf_data, const float3 radiance_out) {
    float3 radiance = radiance_out;
    radiance = radiance + brdf_data->emission;
    radiance = radiance * brdf_data->albedo;
    return radiance;
}

// Gridwalking with 3D-DDA. GPU Gems 3, Chapter 22
// Modified to work on world positions instead of cell positions

// 3D-DDA preparation 
void dda_prep(ray current_ray, float3* step, float3* tmax, float3* tdelta) {
	int originCell = GRID_ID(current_ray.origin);
	//float3 originPos = floor(current_ray.origin / GRID_RES) * GRID_RES + GRID_RES_HALF;

	float3 originPos = GRID_POS(originCell);
	float3 cellMin = originPos - (float3)(GRID_RES_HALF, GRID_RES_HALF, GRID_RES_HALF);  
	float3 cellMax = originPos + (float3)(GRID_RES_HALF, GRID_RES_HALF, GRID_RES_HALF);  
	float3 tmaxNeg = (cellMin - current_ray.origin) / current_ray.direction;  
	float3 tmaxPos = (cellMax - current_ray.origin) / current_ray.direction;  
	*tmax = (current_ray.direction < 0) ? tmaxNeg : tmaxPos;  
	*step = (current_ray.direction < 0) ? (float3)(-1.0f, -1.0f, -1.0f) : (float3)(1.0f, 1.0f, 1.0f);
	*tdelta = fabs((float3)(GRID_RES, GRID_RES, GRID_RES) / current_ray.direction);  
}

// 3D-DDA stepping
void dda_step(float3* pos, float3 step, float3* tmax, float3 tdelta) {
	if((*tmax).x < (*tmax).y) {  
		if((*tmax).x < (*tmax).z) {  
			(*pos).x += step.x;  
			(*tmax).x += tdelta.x;  
		}
		else {  
			(*pos).z += step.z;  
			(*tmax).z += tdelta.z;  
		}  
	}  
	else {  
		if((*tmax).y < (*tmax).z) {  
			(*pos).y += step.y;  
			(*tmax).y += tdelta.y;  
		}  
		else {  
			(*pos).z += step.z;  
			(*tmax).z += tdelta.z;  
		}  
	}
}

// Tracing. Put a ray in, get radiance along ray out. Sample a lot for better results.
inline float3 shade(
	random_state *r, 
	ray ray_in, 
	const __global object* scene, 
	float* ray_probabilities, 
	int* hit_objects,
	const __global int* offsets
) {
    // Iterate over rays until russian roulette decides it's over or we run out of stack
    int depth = 0;
    ray current_ray = ray_in;
	hit best_hit;
	ray_probabilities[0] = ray_in.terminate;
    while(depth < MAX_DEPTH && sample_unit(r) <= current_ray.terminate) {   
        best_hit.hit = 0;
        best_hit.distance = FLT_MAX;
        best_hit.object = 0;
        
		// 3D-DDA variables
		float3 step;
		float3 tmax;
		float3 tdelta;
		dda_prep(current_ray, &step, &tmax, &tdelta);
		
		// Step the ray
		int origin_cell = GRID_ID(current_ray.origin);
		float3 current_pos = (float3)(
			(origin_cell % GRID_RADIUS_2),
			(origin_cell / GRID_RADIUS_2) % GRID_RADIUS_2,
			origin_cell / (GRID_RADIUS_2 * GRID_RADIUS_2)
		);

		// float3 current_pos = current_ray.origin / GRID_RES;
		while(
			fabs(current_pos.x - GRID_RADIUS) < GRID_RADIUS &&
			fabs(current_pos.y - GRID_RADIUS) < GRID_RADIUS &&
			fabs(current_pos.z - GRID_RADIUS) < GRID_RADIUS
		) {
			int currentCell = GRID_CELL(current_pos.x, current_pos.y, current_pos.z);
			int gridPos = offsets[currentCell];
			if(gridPos >= 0) {
				while(scene[gridPos].object_data.n1.x > -5.0f) {
					hit current_hit = triangle_intersect(&scene[gridPos].object_data, current_ray);
					if(current_hit.hit == 1 && current_hit.distance < best_hit.distance) {
						best_hit = current_hit;
						best_hit.object = gridPos;
					}
					gridPos++;
				}

				if(best_hit.hit == 1) {
					if(GRID_CELL(best_hit.position.x, best_hit.position.y, best_hit.position.z) == currentCell) {
						break;
					}
				}
			}
			dda_step(&current_pos, step, &tmax, tdelta);
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

__kernel void ClearImage(
	__write_only image2d_t outImage
) {
	unsigned int i = get_global_id(0);
	if(i < IMAGE_WIDTH * IMAGE_HEIGHT) {
		int coord_x_abs = i % IMAGE_WIDTH;
		int coord_y_abs = floor((float)i / (float)IMAGE_WIDTH);
		write_imagef(outImage, (int2){coord_x_abs, coord_y_abs}, (float4)(0, 0, 0, 0));
	}
}
__kernel void PathTracer(
	__read_only image2d_t inImage, 
	__write_only image2d_t outImage,
	__global int* randInit,
	const __global object* scene,
	const __global int* offsets,
	const float4 cameraPos,
	const float4 cameraUp,
	const float4 cameraFront
) {
	// Coordinates
	unsigned int i = get_global_id(0);
	if(i < IMAGE_WIDTH * IMAGE_HEIGHT) {
		int coord_x_abs = i % IMAGE_WIDTH;
		int coord_y_abs = floor((float)i / (float)IMAGE_WIDTH);
		int coord_x = coord_x_abs - IMAGE_WIDTH / 2;
		int coord_y = -(coord_y_abs - IMAGE_HEIGHT / 2);

		// Randomness
		random_state randstate;
		seed_random(&randstate, randInit[i]);

		// Prepare hit stack
		float ray_probabilities[MAX_DEPTH];
		int hit_objects[MAX_DEPTH];
	
		// Trace
		float3 result = (float3)(0, 0, 0);
		ray trace_ray;
		for(int i = 0; i < SAMPLES; i++) {
			// Generate ray with fuzzed pixel coordinates for AA.
			trace_ray.origin = cameraPos.xyz;
			trace_ray.direction = normalize((float3)(
				(coord_x + sample_unit(&randstate) - 0.5f) * ((float)IMAGE_WIDTH/(float)IMAGE_HEIGHT),
				(coord_y + sample_unit(&randstate) - 0.5f) * ((float)IMAGE_WIDTH/(float)IMAGE_HEIGHT),
				-(float)(IMAGE_WIDTH-IMAGE_WIDTH/2.5f)
			));
			trace_ray.direction =
				trace_ray.direction.x * cross(cameraUp.xyz, cameraFront.xyz) + 
				trace_ray.direction.y * cameraUp.xyz + 
				trace_ray.direction.z * cameraFront.xyz;

			trace_ray.terminate = 1.0f;	
			result += shade(&randstate, trace_ray, scene, ray_probabilities, hit_objects, offsets);
		}
	
		float4 pixelColor = read_imagef(inImage, sampler, (int2){coord_x_abs, coord_y_abs});
		//float4 pixelColor;
		pixelColor += (float4){result.x, result.y, result.z, 1.0f};
		write_imagef(outImage, (int2){coord_x_abs, coord_y_abs}, pixelColor);

		// Write back random
		randInit[i] = randstate.b;
	}
}

//////////////////////////////////////// Filter input calculation //////////////////////////////////////
// Tracing for normals and depth
inline float4 traceNormDepth(
	ray ray_in, 
	const __global object* scene, 
	const __global int* offsets
) {
    ray current_ray = ray_in;
	hit best_hit;
	best_hit.hit = 0;
	best_hit.distance = FLT_MAX;
	best_hit.object = 0;
        
	// 3D-DDA variables
	float3 step;
	float3 tmax;
	float3 tdelta;
	dda_prep(current_ray, &step, &tmax, &tdelta);
		
	// Step the ray
	int origin_cell = GRID_ID(current_ray.origin);
	float3 current_pos = (float3)(
		(origin_cell % GRID_RADIUS_2),
		(origin_cell / GRID_RADIUS_2) % GRID_RADIUS_2,
		origin_cell / (GRID_RADIUS_2 * GRID_RADIUS_2)
	);

	while(
		fabs(current_pos.x - GRID_RADIUS) < GRID_RADIUS &&
		fabs(current_pos.y - GRID_RADIUS) < GRID_RADIUS &&
		fabs(current_pos.z - GRID_RADIUS) < GRID_RADIUS
	) {
		int currentCell = GRID_CELL(current_pos.x, current_pos.y, current_pos.z);
		int gridPos = offsets[currentCell];
		if(gridPos >= 0) {
			while(scene[gridPos].object_data.n1.x > -5.0f) {
				hit current_hit = triangle_intersect(&scene[gridPos].object_data, current_ray);
				if(current_hit.hit == 1 && current_hit.distance < best_hit.distance) {
					best_hit = current_hit;
					best_hit.object = gridPos;
				}
				gridPos++;
			}

			if(best_hit.hit == 1) {
				if(GRID_CELL(best_hit.position.x, best_hit.position.y, best_hit.position.z) == currentCell) {
					break;
				}
			}
		}
		dda_step(&current_pos, step, &tmax, tdelta);
	}

    return((float4)(best_hit.normal.x, best_hit.normal.y, best_hit.normal.z, best_hit.distance));
}


__kernel void NormDepthTracer(
	__write_only image2d_t normDepthImage,
	__global int* randInit,
	const __global object* scene,
	const __global int* offsets,
	const float4 cameraPos,
	const float4 cameraUp,
	const float4 cameraFront
) {
	// Coordinates
	unsigned int i = get_global_id(0);
	if(i < IMAGE_WIDTH * IMAGE_HEIGHT) {
		int coord_x_abs = i % IMAGE_WIDTH;
		int coord_y_abs = floor((float)i / (float)IMAGE_WIDTH);
		int coord_x = coord_x_abs - IMAGE_WIDTH / 2;
		int coord_y = -(coord_y_abs - IMAGE_HEIGHT / 2);

		// Randomness
		random_state randstate;
		seed_random(&randstate, randInit[i]);

		// Trace
		float4 result = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
		ray trace_ray;
		for(int i = 0; i < NORMDEPTH_SAMPLES; i++) {
			// Generate ray with fuzzed pixel coordinates for AA.
			trace_ray.origin = cameraPos.xyz;
			trace_ray.direction = normalize((float3)(
				(coord_x + sample_unit(&randstate) - 0.5f) * ((float)IMAGE_WIDTH/(float)IMAGE_HEIGHT),
				(coord_y + sample_unit(&randstate) - 0.5f) * ((float)IMAGE_WIDTH/(float)IMAGE_HEIGHT),
				-(float)(IMAGE_WIDTH-IMAGE_WIDTH/2.5f)
			));
			trace_ray.direction =
				trace_ray.direction.x * cross(cameraUp.xyz, cameraFront.xyz) + 
				trace_ray.direction.y * cameraUp.xyz + 
				trace_ray.direction.z * cameraFront.xyz;

			trace_ray.terminate = 1.0f;	
			result += traceNormDepth(trace_ray, scene, offsets);
		}
	
		float4 pixelColor = (float4){result.x, result.y, result.z, result.w} / NORMDEPTH_SAMPLES;
		write_imagef(normDepthImage, (int2){coord_x_abs, coord_y_abs}, pixelColor);
	}
}

//////////////////////////////////////// Filtering /////////////////////////////////////////////////////

#define DEPTH_THRESHOLD	0.025f
#define NORM_THRESHOLD	0.9f
#define H_GROUPSIZE_X 32
#define H_GROUPSIZE_Y 4
#define V_GROUPSIZE_X 32
#define V_GROUPSIZE_Y 4
#define KERNEL_RADIUS 4
#define KERNEL_LENGTH (2 * KERNEL_RADIUS + 1)
#define H_RESULT_STEPS 4
#define V_RESULT_STEPS 4

// These functions define discontinuities
bool IsNormalDiscontinuity(float4 n1, float4 n2){
	return fabs(dot(n1, n2)) < NORM_THRESHOLD;
}

bool IsDepthDiscontinuity(float d1, float d2){
	return fabs(d1 - d2) > DEPTH_THRESHOLD;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Horizontal discontinuity filter

//require matching work-group size
__kernel __attribute__((reqd_work_group_size(H_GROUPSIZE_X, H_GROUPSIZE_Y, 1)))
void DiscontinuityHorizontal(
			__global int* d_Disc,
			__read_only image2d_t normalDepthImage
) {
	int Width = IMAGE_WIDTH;
	int Height = IMAGE_HEIGHT;
	int Pitch = IMAGE_WIDTH;

	// We even load unused pixels to the halo area, to keep the code and local memory access simple.
	// Since these loads are coalesced, they introduce no overhead, except for slightly redundant local memory allocation.
	// Each work-item loads H_RESULT_STEPS values + 2 halo values
	// We split the float4 (normal + depth) into an array of float3 and float to avoid bank conflicts.
	__local float tileNormX[H_GROUPSIZE_Y][(H_RESULT_STEPS + 2) * H_GROUPSIZE_X];
	__local float tileNormY[H_GROUPSIZE_Y][(H_RESULT_STEPS + 2) * H_GROUPSIZE_X];
	__local float tileNormZ[H_GROUPSIZE_Y][(H_RESULT_STEPS + 2) * H_GROUPSIZE_X];
	__local float tileDepth[H_GROUPSIZE_Y][(H_RESULT_STEPS + 2) * H_GROUPSIZE_X];

	const int groupSize = H_GROUPSIZE_X * H_RESULT_STEPS;
	const int groupX = get_global_id(0) / H_GROUPSIZE_X;
	const int baseX =  groupX * groupSize;
	const int y = get_global_id(1);
	const int localBaseX = get_local_id(0);
	const int localY = get_local_id(1);
	const int maxX = (H_RESULT_STEPS + 2) * H_GROUPSIZE_X;
	const int maxProcX = H_RESULT_STEPS * H_GROUPSIZE_X;

	// Combined load, seemed unproblematic
	#pragma unroll
	for(int localX = localBaseX; localX < maxX; localX += H_GROUPSIZE_X) {
		int loadX = baseX + localX - H_GROUPSIZE_X;
		if(loadX >= 0 && loadX < Width) {
			float4 nd = read_imagef(normalDepthImage, sampler, (int2){loadX, y});
			tileNormX[localY][localX] = nd.x;
			tileNormY[localY][localX] = nd.y;
			tileNormZ[localY][localX] = nd.z;
			tileDepth[localY][localX] = nd.w;
		}
		else {
			tileNormX[localY][localX] = 0;
			tileNormY[localY][localX] = 0;
			tileNormZ[localY][localX] = 0;
			tileDepth[localY][localX] = 0;
		}
	}

	// Sync threads
	barrier(CLK_LOCAL_MEM_FENCE);

	#pragma unroll
	for(int localX = localBaseX; localX < maxProcX; localX += H_GROUPSIZE_X) {
		if(baseX + localX < Width) {
			int flag = 0;
			float myDepth = tileDepth[localY][localX + H_GROUPSIZE_X];
			float4 myNorm = (float4)(
				tileNormX[localY][localX + H_GROUPSIZE_X],
				tileNormY[localY][localX + H_GROUPSIZE_X],
				tileNormZ[localY][localX + H_GROUPSIZE_X],
				0.0f
			);

			float leftDepth = tileDepth[localY][localX + H_GROUPSIZE_X - 1];
			float4 leftNorm = (float4)(
				tileNormX[localY][localX + H_GROUPSIZE_X - 1],
				tileNormY[localY][localX + H_GROUPSIZE_X - 1],
				tileNormZ[localY][localX + H_GROUPSIZE_X - 1],
				0.0f
			);

			if(IsDepthDiscontinuity(myDepth, leftDepth) || IsNormalDiscontinuity(myNorm, leftNorm)) {
				flag += 1;
			}

			float rightDepth = tileDepth[localY][localX + H_GROUPSIZE_X + 1];
			float4 rightNorm = (float4)(
				tileNormX[localY][localX + H_GROUPSIZE_X + 1],
				tileNormY[localY][localX + H_GROUPSIZE_X + 1],
				tileNormZ[localY][localX + H_GROUPSIZE_X + 1],
				0.0f
			);
			
			if(IsDepthDiscontinuity(myDepth, rightDepth) || IsNormalDiscontinuity(myNorm, rightNorm)) {
				flag += 2;
			}

			d_Disc[y * Pitch + baseX + localX] = flag;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Vertical discontinuity filter

//require matching work-group size
__kernel __attribute__((reqd_work_group_size(V_GROUPSIZE_X, V_GROUPSIZE_Y, 1)))
void DiscontinuityVertical(
			__global int* d_Disc,
			__read_only image2d_t normalDepthImage
) {
	int Width = IMAGE_WIDTH;
	int Height = IMAGE_HEIGHT;
	int Pitch = IMAGE_WIDTH;

	// Local memory
	__local float tileNormX[(V_RESULT_STEPS + 2) * H_GROUPSIZE_Y][H_GROUPSIZE_X];
	__local float tileNormY[(V_RESULT_STEPS + 2) * H_GROUPSIZE_Y][H_GROUPSIZE_X];
	__local float tileNormZ[(V_RESULT_STEPS + 2) * H_GROUPSIZE_Y][H_GROUPSIZE_X];
	__local float tileDepth[(V_RESULT_STEPS + 2) * H_GROUPSIZE_Y][H_GROUPSIZE_X];

	const int groupSize = V_GROUPSIZE_Y * V_RESULT_STEPS;
	const int groupY = get_global_id(1) / V_GROUPSIZE_Y;
	const int baseY =  groupY * groupSize;
	const int x = get_global_id(0);
	const int localBaseY = get_local_id(1);
	const int localX = get_local_id(0);
	const int maxY = (V_RESULT_STEPS + 2) * V_GROUPSIZE_Y;
	const int maxProcY = V_RESULT_STEPS * V_GROUPSIZE_Y;

	// Combined load, seemed unproblematic
	#pragma unroll
	for(int localY = localBaseY; localY < maxY; localY += V_GROUPSIZE_Y) {
		int loadY = baseY + localY - V_GROUPSIZE_Y;
		if(loadY >= 0 && loadY < Height && localY > V_GROUPSIZE_Y - KERNEL_RADIUS - 1 && localY < maxY - V_GROUPSIZE_Y + KERNEL_RADIUS) {
			float4 nd = read_imagef(normalDepthImage, sampler, (int2){x, loadY});
			tileNormX[localY][localX] = nd.x;
			tileNormY[localY][localX] = nd.y;
			tileNormZ[localY][localX] = nd.z;
			tileDepth[localY][localX] = nd.w;
		}
		else {
			tileNormX[localY][localX] = 0;
			tileNormY[localY][localX] = 0;
			tileNormZ[localY][localX] = 0;
			tileDepth[localY][localX] = 0;
		}
	}
	
	// Sync the work-items after loading
	barrier(CLK_LOCAL_MEM_FENCE);

	// Convolve and store the result
	#pragma unroll
	for(int localY = localBaseY; localY < maxProcY; localY += V_GROUPSIZE_Y) {
		if(baseY + localY < Height) {
			int flag = 0;
			float myDepth = tileDepth[localY + V_GROUPSIZE_Y][localX];
			float4 myNorm = (float4)(
				tileNormX[localY + V_GROUPSIZE_Y][localX],
				tileNormY[localY + V_GROUPSIZE_Y][localX],
				tileNormZ[localY + V_GROUPSIZE_Y][localX],
				0.0f
			);

			float upDepth = tileDepth[localY + V_GROUPSIZE_Y - 1][localX];
			float4 upNorm = (float4)(
				tileNormX[localY + V_GROUPSIZE_Y - 1][localX],
				tileNormY[localY + V_GROUPSIZE_Y - 1][localX],
				tileNormZ[localY + V_GROUPSIZE_Y - 1][localX],
				0.0f
			);

			if(IsDepthDiscontinuity(myDepth, upDepth) || IsNormalDiscontinuity(myNorm, upNorm)) {
				flag += 4;
			}

			float downDepth = tileDepth[localY + V_GROUPSIZE_Y + 1][localX];
			float4 downNorm = (float4)(
				tileNormX[localY + V_GROUPSIZE_Y + 1][localX],
				tileNormY[localY + V_GROUPSIZE_Y + 1][localX],
				tileNormZ[localY + V_GROUPSIZE_Y + 1][localX],
				0.0f
			);
			
			if(IsDepthDiscontinuity(myDepth, downDepth) || IsNormalDiscontinuity(myNorm, downNorm)) {
				flag += 8;
			}

			d_Disc[(baseY + localY) * Pitch + x] += flag;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Horizontal convolution filter

//require matching work-group size
__kernel __attribute__((reqd_work_group_size(H_GROUPSIZE_X, H_GROUPSIZE_Y, 1)))
void ConvHorizontal(
			__write_only image2d_t destImage,
			__read_only image2d_t sourceImage,
			__global const int* d_Disc
) {
	int Width = IMAGE_WIDTH;
	int Height = IMAGE_HEIGHT;
	int Pitch = IMAGE_WIDTH;

	// Kernel
	const float c_Kernel[] = {
		0.010284844f,	
		0.0417071f,	
		0.113371652f,	
		0.206576619f,	
		0.252313252f,	
		0.206576619f,	
		0.113371652f,	
		0.0417071f,	
		0.010284844f
	};

	// Local memory
	__local float4 tile[H_GROUPSIZE_Y][(H_RESULT_STEPS + 2) * H_GROUPSIZE_X];
	__local int disc[H_GROUPSIZE_Y][(H_RESULT_STEPS + 2) * H_GROUPSIZE_X];

	const int groupSize = H_GROUPSIZE_X * H_RESULT_STEPS;
	const int groupX = get_global_id(0) / H_GROUPSIZE_X;
	const int baseX =  groupX * groupSize;
	const int y = get_global_id(1);
	const int localBaseX = get_local_id(0);
	const int localY = get_local_id(1);
	const int maxX = (H_RESULT_STEPS + 2) * H_GROUPSIZE_X;
	const int maxProcX = H_RESULT_STEPS * H_GROUPSIZE_X;

	// Combined load
	#pragma unroll
	for(int localX = localBaseX; localX < maxX; localX += H_GROUPSIZE_X) {
		int loadX = baseX + localX - H_GROUPSIZE_X;
		if(loadX >= 0 && loadX < Width) {
			tile[localY][localX] = read_imagef(sourceImage, sampler, (int2){loadX, y});
			disc[localY][localX] = d_Disc[y * Pitch + loadX];
		}
		else {
			tile[localY][localX] = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
			disc[localY][localX] = 0;
		}
	}
	
	// Sync the work-items after loading
	barrier(CLK_LOCAL_MEM_FENCE);

	// Convolve and store the result
	#pragma unroll
	for(int localX = localBaseX; localX < maxProcX; localX += H_GROUPSIZE_X) {
		if(baseX + localX < Width) {
			float4 temp = c_Kernel[KERNEL_RADIUS] * tile[localY][localX + H_GROUPSIZE_X];
			float weight = c_Kernel[KERNEL_RADIUS];

			// To the left
			if((disc[localY][localX + H_GROUPSIZE_X] & 1) == 0) {
				for(int i = -1; i >= -KERNEL_RADIUS; i--) {
					temp += c_Kernel[KERNEL_RADIUS + i] * tile[localY][localX + H_GROUPSIZE_X + i];
					weight += c_Kernel[KERNEL_RADIUS + i];
					if((disc[localY][localX + H_GROUPSIZE_X + i] & 1) != 0) {
						break;
					}
				}
			}

			// To the right
			if((disc[localY][localX + H_GROUPSIZE_X] & 2) == 0) {
				for(int i = 1; i <= KERNEL_RADIUS; i++) {
					temp += c_Kernel[KERNEL_RADIUS + i] * tile[localY][localX + H_GROUPSIZE_X + i];
					weight += c_Kernel[KERNEL_RADIUS + i];
					if((disc[localY][localX + H_GROUPSIZE_X + i] & 2) != 0) {
						break;
					}
				}
			}

			// Write back
			write_imagef(destImage, (int2){baseX + localX, y},  temp / weight);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Vertical convolution filter

//require matching work-group size
__kernel __attribute__((reqd_work_group_size(V_GROUPSIZE_X, V_GROUPSIZE_Y, 1)))
void ConvVertical(
			__write_only image2d_t destImage,
			__read_only image2d_t sourceImage,
			__global const int* d_Disc
) {
	int Width = IMAGE_WIDTH;
	int Height = IMAGE_HEIGHT;
	int Pitch = IMAGE_WIDTH;

	// Kernel
	const float c_Kernel[9] = {
		0.010284844f,	
		0.0417071f,	
		0.113371652f,	
		0.206576619f,	
		0.252313252f,	
		0.206576619f,	
		0.113371652f,	
		0.0417071f,	
		0.010284844f
	};

	// Passthrough
	/*for(int i = 0; i < V_RESULT_STEPS; i++) {
		d_Dst[(get_global_id(1) * V_RESULT_STEPS + i) * Pitch + get_global_id(0)] = d_Src[(get_global_id(1) * V_RESULT_STEPS + i) * Pitch + get_global_id(0)];
	}*/

	// Local memory
	__local float4 tile[(V_RESULT_STEPS + 2) * H_GROUPSIZE_Y][H_GROUPSIZE_X];
	__local int disc[(V_RESULT_STEPS + 2) * H_GROUPSIZE_Y][H_GROUPSIZE_X];

	const int groupSize = V_GROUPSIZE_Y * V_RESULT_STEPS;
	const int groupY = get_global_id(1) / V_GROUPSIZE_Y;
	const int baseY =  groupY * groupSize;
	const int x = get_global_id(0);
	const int localBaseY = get_local_id(1);
	const int localX = get_local_id(0);
	const int maxY = (V_RESULT_STEPS + 2) * V_GROUPSIZE_Y;
	const int maxProcY = V_RESULT_STEPS * V_GROUPSIZE_Y;

	// Combined load, seemed unproblematic
	#pragma unroll
	for(int localY = localBaseY; localY < maxY; localY += V_GROUPSIZE_Y) {
		int loadY = baseY + localY - V_GROUPSIZE_Y;
		if(loadY >= 0 && loadY < Height && localY > V_GROUPSIZE_Y - KERNEL_RADIUS - 1 && localY < maxY - V_GROUPSIZE_Y + KERNEL_RADIUS) {
			tile[localY][localX] = read_imagef(sourceImage, sampler, (int2){x, loadY});
			disc[localY][localX] = d_Disc[loadY * Pitch + x];
		}
		else {
			tile[localY][localX] = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
			disc[localY][localX] = 0;
		}
	}
	
	// Sync the work-items after loading
	barrier(CLK_LOCAL_MEM_FENCE);

	// Convolve and store the result
	#pragma unroll
	for(int localY = localBaseY; localY < maxProcY; localY += V_GROUPSIZE_Y) {
		if(baseY + localY < Height) {
			float4 temp = c_Kernel[KERNEL_RADIUS] * tile[localY + V_GROUPSIZE_Y][localX];
			float weight = c_Kernel[KERNEL_RADIUS];

			// Gunbuster
			if((disc[localY + V_GROUPSIZE_Y][localX] & 4) == 0) {
				for(int i = -1; i >= -KERNEL_RADIUS; i--) {
					temp += c_Kernel[KERNEL_RADIUS + i] * tile[localY + V_GROUPSIZE_Y + i][localX];
					weight += c_Kernel[KERNEL_RADIUS + i];
					if((disc[localY + V_GROUPSIZE_Y + i][localX] & 4) != 0) {
						break;
					}
				}
			}

			// To the bottom
			if((disc[localY + V_GROUPSIZE_Y][localX] & 8) == 0) {
				for(int i = 1; i <= KERNEL_RADIUS; i++) {
					temp += c_Kernel[KERNEL_RADIUS + i] * tile[localY + V_GROUPSIZE_Y + i][localX];
					weight += c_Kernel[KERNEL_RADIUS + i];
					if((disc[localY + V_GROUPSIZE_Y + i][localX] & 8) != 0) {
						break;
					}
				}
			}

			// Write back
			write_imagef(destImage, (int2){x, (baseY + localY)}, temp / weight);
		}
	}
}
