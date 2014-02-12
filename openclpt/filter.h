#ifndef __FILTER_H__
#define __FILTER_H__

#include "clgl.h"

void prepareFilter(cl_program program, int width, int height);
GLuint prepareNormalDepth(cl_mem randInit, cl_mem scene, cl_mem offsets, float* pos, float* up, float* front);
GLuint runFilterKernels(cl_mem inImage);

#endif