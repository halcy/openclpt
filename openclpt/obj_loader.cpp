/**
 * Loads a subset and variant of wavefront obj
 *
 * The subset being triangles and per-vertex normals, and the variant
 * being material names must be numbers corresponding to the materials
 * file.
 * 
 * There is no error handling whatsoever, so please stick to the format.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_structures.h"

#define MAX_VERTICES 400000
#define MAX_TRIANGLES 400000
#define MAX_LINE_LEN 1024
#define MAX_MATERIALS 100

// "THIS FUNCTION IS UNSAFE"
#pragma warning(disable: 4996)

vec3_t readFloatVector(char* source) {
	vec3_t v;
	int ret = sscanf(source, "%f %f %f", &v.x, &v.y, &v.z);
	//printf("Loaded (%d) vector: %f %f %f\n", ret, v.x, v.y, v.z);
	v.w = 1.0f;
	return v;
}

simple_brdf_data* loadMaterials(char* fileName) {
	// Alloc memory for materials
	simple_brdf_data* materials = (simple_brdf_data*)malloc(sizeof(simple_brdf_data) * MAX_MATERIALS);
	
	// Read materials
	FILE* objFile = fopen (fileName, "rt");
	char line[MAX_LINE_LEN];
	int currentMaterial = 0;

	while(fgets(line, MAX_LINE_LEN, objFile) != NULL) {
		if(line[0] == 'm') {
			sscanf(&line[2], "%d", &currentMaterial);
		}

		if(line[0] == 'e') {
			materials[currentMaterial].emission = readFloatVector(&line[2]);
		}

		if(line[0] == 'a') {
			materials[currentMaterial].albedo = readFloatVector(&line[2]);
		}

		if(line[0] == 'p') {
			vec3_t properties = readFloatVector(&line[2]);
			materials[currentMaterial].reflectivity = properties.x;
			materials[currentMaterial].transparency = properties.y;
			materials[currentMaterial].specularity = properties.z;
		}
	}

	return materials;
}

object* loadObj(char* fileName, int* numTriangles, char* materialsFileName) {
	// First, load materials
	simple_brdf_data* materials = loadMaterials(materialsFileName);

	// Set up data for loading
	vec3_t* vertices = (vec3_t*)malloc(sizeof(vec3_t) * MAX_VERTICES);
	vec3_t* normals = (vec3_t*)malloc(sizeof(vec3_t) * MAX_VERTICES);
	int* vertexIndices = (int*)malloc(sizeof(int) * MAX_TRIANGLES * 3);
	int* normalIndices = (int*)malloc(sizeof(int) * MAX_TRIANGLES * 3);
	int* materialIndices = (int*)malloc(sizeof(int) * MAX_TRIANGLES);

	int vertCount = 0;
	int normCount = 0;
	int triCount = 0;
	int currentMaterial = 0;

	// Load
	FILE* objFile = fopen (fileName, "rt");
	char line[MAX_LINE_LEN];
	while(fgets(line, MAX_LINE_LEN, objFile) != NULL) {
		if(line[0] == 'v' && line[1] != 'n') {
			vertices[vertCount] = readFloatVector(&line[2]);
			vertices[vertCount].x = -vertices[vertCount].x;
			vertCount++;
		}

		if(line[0] == 'v' && line[1] == 'n') {
			normals[normCount] = readFloatVector(&line[3]);
			normals[normCount].x = -normals[normCount].x;
			normCount++;
		}

		if(line[0] == 'f') {
			int v1, v2, v3, n1, n2, n3;
			//printf("%s\n", line);
			sscanf(&line[2], "%d//%d %d//%d %d//%d", &v1, &n1, &v2, &n2, &v3, &n3);
			//printf("Loaded indices: %d//%d %d//%d %d//%d\n", v1, n1, v2, n2, v3, n3);
			vertexIndices[triCount * 3] = v1 - 1;
			vertexIndices[triCount * 3 + 1] = v2 - 1;
			vertexIndices[triCount * 3 + 2] = v3 - 1;
			normalIndices[triCount * 3] = n1 - 1;
			normalIndices[triCount * 3 + 1] = n2 - 1;
			normalIndices[triCount * 3 + 2] = n3 - 1;
			materialIndices[triCount] = currentMaterial;
			triCount++;
		}

		if(line[0] == 'u') {
			sscanf(&line[7], "%d", &currentMaterial);
		}
	}
	printf("%s: Loaded %d vertices, %d normals, %d faces.\n", fileName, vertCount, normCount, triCount);

	// Convert to internal format (note: rewinding tris to ccw)
	object* scene = (object*)malloc(sizeof(object) * triCount);
	for(int i = 0; i < triCount; i++) {
		scene[i].object_data.p1 = vertices[vertexIndices[i * 3]];
		scene[i].object_data.p2 = vertices[vertexIndices[i * 3 + 2]];
		scene[i].object_data.p3 = vertices[vertexIndices[i * 3 + 1]];
		scene[i].object_data.n1 = normals[normalIndices[i * 3]];
		scene[i].object_data.n2 = normals[normalIndices[i * 3 + 2]];
		scene[i].object_data.n3 = normals[normalIndices[i * 3 + 1]];
		scene[i].brdf_data = materials[materialIndices[i]];
		/*printf(
			"Tri %d: %f %f %f / %f %f %f\n", 
			i,
			scene[i].object_data.p1.x,
			scene[i].object_data.p1.z,
			scene[i].object_data.p1.z,
			scene[i].object_data.n1.x,
			scene[i].object_data.n1.y,
			scene[i].object_data.n1.z
		);*/
	}

	// Free and return
	free(vertices);
	free(normals);
	free(vertexIndices);
	free(normalIndices);
	free(materialIndices);
	free(materials);

	printf("Done\n");

	*numTriangles = triCount;
	return scene;
}
