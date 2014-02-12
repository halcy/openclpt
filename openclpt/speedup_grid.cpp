#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define GRID_ELEMENTS_EXPECTED_MULT 4.0f
#define GRID_ELEMENTS_MIN 800000

#include "data_structures.h"

void maxProjection(vec3_t* points, int count, vec3_t axis, float* min, float* max) {
	*min = 1000000000.0f;
	*max = -1000000000.0f;
	// printf("Axis: %f %f %f\n", axis.x, axis.y, axis.z);
	for(int i = 0; i < count; i++) {
		float val = dot(axis, points[i]);
		// printf("Dot is %f\n", val);
		if(val < *min) {
			// printf("new min\n");
			*min = val;
		}
		if(val > *max) {
			// printf("new max\n");
			*max = val;
		}
	}
}

// Ported from http://stackoverflow.com/questions/17458562/efficient-aabb-triangle-intersection-in-c-sharp
bool triangleIntersectsAABB(aabb_data box, triangle_object_data triangle) {
    float triangleMin;
	float triangleMax;
    float boxMin;
	float boxMax;

	// Find triangle face normal
	vec3_t triangleNorm = norm(cross(sub(triangle.p2, triangle.p1), sub(triangle.p3, triangle.p1)));

	// Construct vertex arrays
	vec3_t triangleVerts[3];
	triangleVerts[0] = triangle.p1;
	triangleVerts[1] = triangle.p2;
	triangleVerts[2] = triangle.p3;

    // Test the box normals (x-, y- and z-axes)
	vec3_t boxNormals[3] = {
		vec3(1.0f, 0.0f, 0.0f),
		vec3(0.0f, 1.0f, 0.0f),
		vec3(0.0f, 0.0f, 1.0f)
	};

    for(int i = 0; i < 3; i++) {
        vec3_t n = boxNormals[i];
        maxProjection(triangleVerts, 3, boxNormals[i], &triangleMin, &triangleMax);
		boxMin = i == 0 ? box.vertices[0].x : (i == 1 ? box.vertices[0].y : box.vertices[0].z);
		boxMax = i == 0 ? box.vertices[7].x : (i == 1 ? box.vertices[7].y : box.vertices[7].z);
		if (triangleMax < boxMin || triangleMin > boxMax) {
			//printf("Boxnorm hit\n");
            return false;
		}
    }

    // Test the triangle normal
    float triangleOffset = dot(triangleNorm, triangle.p1);
	maxProjection(box.vertices, 8, triangleNorm, &boxMin, &boxMax);
    if (boxMax < triangleOffset || boxMin > triangleOffset) {
		// printf("Trinorm hit\n");
        return false; // No intersection possible.
	}

    // Test the nine edge cross-products
	vec3_t triangleEdges[3];
	triangleEdges[0] = sub(triangle.p1, triangle.p2);
	triangleEdges[1] = sub(triangle.p2, triangle.p3);
	triangleEdges[2] = sub(triangle.p3, triangle.p1);
	
    for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			// The box normals are the same as it's edge tangents
			vec3_t axis = cross(triangleEdges[i], boxNormals[j]);
			if(length(axis) > 0.000001f) {
				maxProjection(box.vertices, 8, axis, &boxMin, &boxMax);
				maxProjection(triangleVerts, 3, axis, &triangleMin, &triangleMax);
				if (boxMax < triangleMin || boxMin > triangleMax) {
					// printf("Edge crossprod hit %d %d -> %f %f %f %f\n", i, j, boxMax, triangleMin, boxMin, triangleMax);	
					return false; // No intersection possible
				}
			}
		}
	}

    // No separating axis found.
    return true;
}

// Find maximum free cells starting from some cell
int cell_diff[26][3] = {
	{-1, -1, -1},
	{-1, -1,  0},
	{-1, -1,  1},
	{-1,  0, -1},
	{-1,  0,  0},
	{-1,  0,  1},
	{-1,  1, -1},
	{-1,  1,  0},
	{-1,  1,  1},
	
	{ 0, -1, -1},
	{ 0, -1,  0},
	{ 0, -1,  1},
	{ 0,  0, -1},
	{ 0,  0,  1},
	{ 0,  1, -1},
	{ 0,  1,  0},
	{ 0,  1,  1},
	
	{ 1, -1, -1},
	{ 1, -1,  0},
	{ 1, -1,  1},
	{ 1,  0, -1},
	{ 1,  0,  0},
	{ 1,  0,  1},
	{ 1,  1, -1},
	{ 1,  1,  0},
	{-1,  1,  1}
};

int free_cells(int* offset, bool* mark, int x, int y, int z, int extent2) {
	// BFS queue
	int testCells[2048][4];
	int testPointer = 0;
	
	// Prepare BFS
	mark[x + y * extent2 + z * extent2 * extent2] = true;
	testCells[testPointer][0] = x;
	testCells[testPointer][1] = y;
	testCells[testPointer][2] = z;
	testCells[testPointer][3] = 1;
	testPointer++;
	testPointer = testPointer % 2048;

	// Perform BFS
	//printf("Go BFS\n");
	for(int i = 0; i < testPointer; i++) {
		
		int currentDepth = testCells[i][3];
		int x1 = testCells[i][0];
		int y1 = testCells[i][1];
		int z1 = testCells[i][2];
		//printf("%d -> %d %d %d -> %d\n", i, x1, y1, z1, currentDepth);

		// Bounds check
		if(x1 < 0 || x1 >= extent2) return currentDepth;
		if(y1 < 0 || y1 >= extent2) return currentDepth;
		if(z1 < 0 || z1 >= extent2) return currentDepth;

		// Fill check
		int id1 = x1 + y1 * extent2 + z1 * extent2 * extent2;
		if(offset[id1] >= 0) return currentDepth;

		// Cell is free, mark and add unmarked neighbours to queue
		for(int j = 0; j < 26; j++) {
			int x2 = x1 + cell_diff[j][0];
			int y2 = y1 + cell_diff[j][1];
			int z2 = z1 + cell_diff[j][2];
			if(x2 < 0 && x2 >= extent2 && y < 0 && y >= extent2 && z < 0 && z >= extent2) {
				int id2 = x2 + y2 * extent2 + z2 * extent2 * extent2;
				if(mark[id2] == false) {
					mark[id2] = true;
					testCells[testPointer][0] = x2;
					testCells[testPointer][1] = y2;
					testCells[testPointer][2] = z2;
					testCells[testPointer][3] = currentDepth + 1;
					testPointer++;
					testPointer = testPointer % 2048;
				}
			}
		}
	}

	// This should never actually happen - be conservative if it does
	return 1;
}

// Ideally this would be OpenCL, too.
// Also ideally this would have some form of deduplication
int createGrid(object* scene, int objects, int extent, float res, bool is_sub, object** gridSceneRet, int** offsetsRet) {
	// End marker
	object endMarker;
	endMarker.object_data.n1 = vec3(-10.0f, 0.0f, 0.0f); // Super dangerous

	// Recurse for lower level grid
	object* lowerGrid;
	int* lowerOffsets;
	int maxElements;
	if(extent > 1) {
		maxElements = GRID_ELEMENTS_EXPECTED_MULT * createGrid(scene, objects, extent / 2, res * 2.0f, true, &lowerGrid, &lowerOffsets);
	}
	else {
		maxElements = GRID_ELEMENTS_EXPECTED_MULT * objects;
	}
	maxElements = maxElements < GRID_ELEMENTS_MIN ? GRID_ELEMENTS_MIN : maxElements;

	// Pre-allocate grid memory
	int extent2 = 2 * extent;
	printf("Allocating max: %d for extent level %d\n", maxElements, extent);
	object* gridScene = (object*)malloc((int)(sizeof(object) * maxElements));
	int* offsets = (int*)malloc(sizeof(int) * extent2 * extent2 * extent2);

	// Do griding
	int currentIndex = 0;
	for(int x = 0; x < extent2; x++) {
		printf("Grid extent %d x layer %d computing...\n", extent, x);
		for(int y = 0; y < extent2; y++) {
			for(int z = 0; z < extent2; z++) {
				// Save start index for later
				int startIndex = currentIndex;

				// Build box
				vec3_t center = vec3(
					((float)(x - extent) + 0.5f) * res,
					((float)(y - extent) + 0.5f) * res,
					((float)(z - extent) + 0.5f) * res
				);
				vec3_t radius = vec3(res * 0.5f, res * 0.5f, res * 0.5f);
				aabb_data box = aabb(center, radius);
				
				int triCount = 0;
				if(extent == 1) {
					for(int i = 0; i < objects; i++) {
						//printf("Call tri %d\n", i);
						if(triangleIntersectsAABB(box, scene[i].object_data)) {
							//printf("Isect write %d\n", currentIndex);
							gridScene[currentIndex] = scene[i];
							currentIndex++;
							triCount++;
						}
					}
				}
				else {
					int offset = lowerOffsets[x/2 + (y/2) * (extent2/2) + (z/2) * (extent2/2) * (extent2/2)];
					//printf("Offset is: %d\n", offset);
					int sa = offset;
					while(offset >= 0 && lowerGrid[offset].object_data.n1.x > -5.0f) {
						if(triangleIntersectsAABB(box, lowerGrid[offset].object_data)) {
							//printf("Isect write %d\n", currentIndex);
							gridScene[currentIndex] = lowerGrid[offset];
							currentIndex++;
							triCount++;
						}
						offset++;
					}
					//printf("Stopped at %d, ran through %d cells\n", offset, offset - sa);
				}

				// Figure out index and store offset
				int boxIndex = x + y * extent2 + z * extent2 * extent2;
				if(triCount > 0) {
					// End marker
					gridScene[currentIndex] = endMarker;
					currentIndex++;

					offsets[boxIndex] = startIndex;
					//printf("%d %d %d [%d]: %d tris @ %d\n", x, y, z, boxIndex, triCount, startIndex);
				}
				else {
					offsets[boxIndex] = -1;
					//printf("%d %d %d [%d]: empty\n", x, y, z, boxIndex);
				}
			}
		}
	}

	/*if(!is_sub) {
		for(int x = 1; x < extent2 - 1; x++) {
			printf("Grid extent %d x layer %d computing empy space skip...\n", extent, x);
			for(int y = 1; y < extent2 - 1; y++) {
				printf("Grid extent %d y layer %d computing empy space skip...\n", extent, y);
				for(int z = 1; z < extent2 - 1; z++) {
					int boxIndex = x + y * extent2 + z * extent2 * extent2;
					if(offsets[boxIndex] <= 0) {
						bool* mark = (bool*)malloc(sizeof(bool) * extent2 * extent2 * extent2);
						for(int i = 0; i < extent2 * extent2 * extent2; i++) {
							mark[i] = false;
						}
						offsets[boxIndex] = -free_cells(offsets, mark, x, y, z, extent2);
						if(offsets[boxIndex] != -1) {
							printf("%d %d %d: %d free\n", x, y, z, offsets[boxIndex]);
						}
						free(mark);
					}
				}
			}
		}
	}*/

	*offsetsRet = offsets;
	*gridSceneRet = gridScene;

	if(extent > 1) {
		free(lowerGrid);
		free(lowerOffsets);
	}

	printf("Triangles for extent: %d\n", currentIndex);
	return currentIndex;
}