#ifndef __DATA_H__
#define __DATA_H__

#include "data_structures.h"
#include "clgl.h"

// Buffer objects for OpenGL and OpenCL
struct {
	GLuint vertexBuffer;
	GLuint elementBuffer;
} screenQuad;

// The display shader
struct {
	GLuint shaderProgram;
	GLuint vertexPosition;

	GLuint textureLocation;
	GLuint iterCount;
	GLuint userScale;

	GLuint texture[2];
} displayShader;

// A bunch of kernels for the tracer
struct {
	cl_program program;

	cl_kernel renderKernel;
	cl_kernel clearKernel;

	int targetBufferIndex;
	float sampleCount;
	cl_mem renderBuffer[2];
} tracer;

// A scene
struct {
	object* scene;
	int triCount;

	object* gridScene;
	int gridSceneSize;
	int* gridOffsets;
	int gridRadius;
	float gridRes;

	bool onDevice;
	cl_mem sceneBuffer;
	cl_mem offsetBuffer;
} scene;

// Camera settings
#define CAM_MOVESPEED 0.02f
#define CAM_ROTSPEED 0.0005f
struct {
	Vector pos;
	Vector front;
	Vector up;
	float elevation;
} camera;

#endif