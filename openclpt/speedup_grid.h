#ifndef __SPEEDUP_GRID_H__
#define __SPEEDUP_GRID_H__

#include "data_structures.h"

bool triangleIntersectsAABB(aabb_data box, triangle_object_data triangle);
int createGrid(object* scene, int objects, int extent, float res, bool is_sub, object** gridSceneRet, int** offsetsRet);

#endif