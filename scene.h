/**
 * Scene: Cornell box with cooler colours and a bunch of spheres.
 *
 * (c) L. Diener, 2013
 */
 
#ifndef __SCENE_H__
#define __SCENE_H__

// How many objects in this scene
#define SCENE_OBJECT_COUNT 13
 
#include "data_structures.h"
 
// (THAT scene.)
inline static void prepare_scene(object* scene, float animation_param) {    
	// Parameters of sorts
    float at_depth = -4.9f;
    float brightness = 500.0f - fmax(0.0f, (-cos(animation_param)+0.3f))*(480.0/1.3f);
	
    // Temporary objects
	sphere_object_data spheres[13];
	simple_brdf_data materials[8];

    // Materials.
    // Bright (warm-)white light
    materials[0].emission = vec3(brightness, brightness*0.85f, brightness*0.8f);
    materials[0].albedo = vec3(1.0f, 1.0f, 1.0f);
    materials[0].reflectivity = 0;
    materials[0].specularity = 0;

    // White diffuse
    materials[1].emission = vec3(0, 0, 0);
    materials[1].albedo = vec3(0.99999f, 0.99999f, 0.99999f);
    materials[1].reflectivity = 0;
    materials[1].specularity = 0;

    // Purple diffuse
    materials[2].emission = vec3(0, 0, 0);
    materials[2].albedo = vec3(0.44444f, 0, 0.55555f);
    materials[2].reflectivity = 0;
    materials[2].specularity = 0;

    // Yellow diffuse
    materials[3].emission = vec3(0, 0, 0);
    materials[3].albedo = vec3(0.99999f, 0.99999f, 0);
    materials[3].reflectivity = 0;
    materials[3].specularity = 0;

    // Imperfect specular reflector
    materials[4].emission = vec3(0, 0, 0);
    materials[4].albedo = vec3(0.99999f, 0.99999f, 0.99999f);
    materials[4].reflectivity = 0.75f;
    materials[4].specularity = 2.5f;

    // Perfect specular reflector
    materials[5].emission = vec3(0, 0, 0);
    materials[5].albedo = vec3(0.99999f, 0.99999f, 0.99999f);
    materials[5].reflectivity = 1.0f;
    materials[5].specularity = 100000.0f;

    // Pretty good coloured specular reflector
    materials[6].emission = vec3(0, 0, 0);
    materials[6].albedo = vec3(0.99999f, 0.79999f, 0.29999f);
    materials[6].reflectivity = 0.6f;
    materials[6].specularity = 15.5f;

    // Coloured light
    materials[7].emission = vec3(0.0f, 1700.0f/3.0f, 1800.0f/3.0f);
    materials[7].albedo = vec3(0.99999f, 0.99999f, 0.99999f);
    materials[7].reflectivity = 0.0f;
    materials[7].specularity = 0.0f;
    
    // Geometry.
    // Light
    spheres[0].position = vec3(0, -5.0f, at_depth);
    spheres[0].radius = 3.3f;
    scene[0].object_data = spheres[0];
    scene[0].brdf_data = materials[0];
    
    // Ceiling
    spheres[1].position = vec3(0, -502.0f, at_depth);
    spheres[1].radius = 500.0f;
    scene[1].object_data = spheres[1];
    scene[1].brdf_data = materials[1];
    
    // Floor
    spheres[2].position = vec3(0, 502.0f, at_depth);
    spheres[2].radius = 500.0f;
    scene[2].object_data = spheres[2];
    scene[2].brdf_data = materials[1];
    
    // Left
    spheres[3].position = vec3(-502.0f, 0, at_depth);
    spheres[3].radius = 500.0f;
    scene[3].object_data = spheres[3];
    scene[3].brdf_data = materials[3];
    
    // Right
    spheres[4].position = vec3(502.0f, 0, at_depth);
    spheres[4].radius = 500.0f;
    scene[4].object_data = spheres[4];
    scene[4].brdf_data = materials[2];
    
    // Back
    spheres[5].position = vec3(0, 0, -502.0f + at_depth);
    spheres[5].radius = 500.0f;
    scene[5].object_data = spheres[5];
    scene[5].brdf_data = materials[1];
    
    // Front - a little farther away because camera.
    spheres[6].position = vec3(0, 0, 506.0f + at_depth);
    spheres[6].radius = 500.0f;
    scene[6].object_data = spheres[6];
    scene[6].brdf_data = materials[1];
    
    // Gratuitous decorative floor spheres
    spheres[7].position = vec3(-0.8f, 1.3f, -0.1f + at_depth);
    spheres[7].radius = 0.7f;
    scene[7].object_data = spheres[7];
    scene[7].brdf_data = materials[1];
    
    spheres[8].position = vec3(1.3f, 1.6f, 0.2f + at_depth);
    spheres[8].radius = 0.4f;
    scene[8].object_data = spheres[8];
    scene[8].brdf_data = materials[5];

    spheres[9].position = vec3(0.4f, 1.45f, -0.9f + at_depth);
    spheres[9].radius = 0.55f;
    scene[9].object_data = spheres[9];
    scene[9].brdf_data = materials[4];

    spheres[10].position = vec3(-0.08f, 1.7f, 0.6f + at_depth);
    spheres[10].radius = 0.3f;
    scene[10].object_data = spheres[10];
    scene[10].brdf_data = materials[6];

	float sphere_shift = (-cos(animation_param)+1.0f)*1.5;
    spheres[11].position = vec3(0.55f, 1.8f - sphere_shift, -0.05f + at_depth);
    spheres[11].radius = 0.2f;
    scene[11].object_data = spheres[11];
    scene[11].brdf_data = materials[7];
	
	// Light leak prevention
	spheres[12].position = vec3(0, -52.0f, at_depth);
    spheres[12].radius = 50.000f;
    scene[12].object_data = spheres[12];
    scene[12].brdf_data = materials[1];
}

#endif