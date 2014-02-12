#include "glhelpers.h"
#include "clgl.h"

// A bunch of kernels and data for the filter
struct {
	cl_program program;
	
	GLuint normalsDepthTexture;
	GLuint tempTexture;
	GLuint resultTexture;

	cl_mem normalsDepthBuffer;
	cl_mem discBuffer;
	cl_mem tempBuffer;
	cl_mem resultBuffer;

	cl_kernel normalDepthKernel;
	cl_kernel discHorizontalKernel;
	cl_kernel discVerticalKernel;
	cl_kernel convHorizontalKernel;
	cl_kernel convVerticalKernel;

	int width;
	int height;
} filter;

void prepareFilter(cl_program program, int width, int height) {
	// Create kernels
	filter.normalDepthKernel = clCreateKernel(program, "NormDepthTracer", NULL);
	filter.discHorizontalKernel = clCreateKernel(program, "DiscontinuityHorizontal", NULL);
	filter.discVerticalKernel = clCreateKernel(program, "DiscontinuityVertical", NULL);
	filter.convHorizontalKernel = clCreateKernel(program, "ConvHorizontal", NULL);
	filter.convVerticalKernel = clCreateKernel(program, "ConvVertical", NULL);

	// Create buffers
	filter.normalsDepthTexture = makeTextureBuffer(width, height, GL_RGBA, GL_RGBA32F);
	filter.normalsDepthBuffer = clCreateFromGLTexture2D(clContext(), CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, filter.normalsDepthTexture, 0);
	filter.tempTexture = makeTextureBuffer(width, height, GL_RGBA, GL_RGBA32F);
	filter.tempBuffer = clCreateFromGLTexture2D(clContext(), CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, filter.tempTexture, 0);
	filter.resultTexture = makeTextureBuffer(width, height, GL_RGBA, GL_RGBA32F);
	filter.resultBuffer = clCreateFromGLTexture2D(clContext(), CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, filter.resultTexture, 0);
	filter.discBuffer = clCreateBuffer(
		clContext(), 
		CL_MEM_READ_WRITE, 
		sizeof(int) * width * height, 
		0, 
		NULL
	);

	// Bind static kernel attributes
	clSetKernelArg(filter.normalDepthKernel, 0, sizeof(cl_mem), (void*)&(filter.normalsDepthBuffer));

	clSetKernelArg(filter.discHorizontalKernel, 0, sizeof(cl_mem), (void*)&(filter.discBuffer));
	clSetKernelArg(filter.discHorizontalKernel, 1, sizeof(cl_mem), (void*)&(filter.normalsDepthBuffer));

	clSetKernelArg(filter.discVerticalKernel, 0, sizeof(cl_mem), (void*)&(filter.discBuffer));
	clSetKernelArg(filter.discVerticalKernel, 1, sizeof(cl_mem), (void*)&(filter.normalsDepthBuffer));

	clSetKernelArg(filter.convHorizontalKernel, 0, sizeof(cl_mem), (void*)&(filter.tempBuffer));
	clSetKernelArg(filter.convHorizontalKernel, 2, sizeof(cl_mem), (void*)&(filter.discBuffer));
		
	clSetKernelArg(filter.convVerticalKernel, 0, sizeof(cl_mem), (void*)&(filter.resultBuffer));
	clSetKernelArg(filter.convVerticalKernel, 1, sizeof(cl_mem), (void*)&(filter.tempBuffer));
	clSetKernelArg(filter.convVerticalKernel, 2, sizeof(cl_mem), (void*)&(filter.discBuffer));

	// Save dimensions
	filter.width = width;
	filter.height = height;
}

GLuint prepareNormalDepth(cl_mem randInit, cl_mem scene, cl_mem offsets, float* pos, float* up, float* front) {
	acquireGLBuffer(filter.normalsDepthBuffer);

	cl_int numPixels = filter.width * filter.height;
	cl_uint workSize[3] = {numPixels, 0, 0};
	cl_uint workgroupSize[3] = {256, 0, 0};

	clSetKernelArg(filter.normalDepthKernel, 1, sizeof(cl_mem), (void*)&randInit);
	clSetKernelArg(filter.normalDepthKernel, 2, sizeof(cl_mem), (void*)&scene);
	clSetKernelArg(filter.normalDepthKernel, 3, sizeof(cl_mem), (void*)&offsets);
	clSetKernelArg(filter.normalDepthKernel, 4, sizeof(cl_float4), (void*)pos);
	clSetKernelArg(filter.normalDepthKernel, 5, sizeof(cl_float4), (void*)up);
	clSetKernelArg(filter.normalDepthKernel, 6, sizeof(cl_float4), (void*)front);
	clRunKernel(filter.normalDepthKernel, workSize, workgroupSize);

	size_t filterLocalSize[2] = {32, 4};
	size_t globalWorkSizeH[2] = {filter.width / 4, filter.height};	
	size_t globalWorkSizeV[2] = {filter.width, filter.height / 4};

	clEnqueueNDRangeKernel(clCommandQueue(), filter.discHorizontalKernel, 2, NULL, globalWorkSizeH, filterLocalSize, 0, NULL, NULL);
	clFinish(clCommandQueue());

	clEnqueueNDRangeKernel(clCommandQueue(), filter.discVerticalKernel, 2, NULL, globalWorkSizeV, filterLocalSize, 0, NULL, NULL);
	clFinish(clCommandQueue());

	releaseGLBuffer(filter.normalsDepthBuffer);

	return filter.normalsDepthTexture;
}

GLuint runFilterKernels(cl_mem inImage) {
	acquireGLBuffer(inImage);
	acquireGLBuffer(filter.normalsDepthBuffer);
	acquireGLBuffer(filter.tempBuffer);
	acquireGLBuffer(filter.resultBuffer);

	size_t filterLocalSize[2] = {32, 4};
	size_t globalWorkSizeH[2] = {filter.width / 4, filter.height};	
	size_t globalWorkSizeV[2] = {filter.width, filter.height / 4};

	clSetKernelArg(filter.convHorizontalKernel, 1, sizeof(cl_mem), (void*)&inImage);
	clEnqueueNDRangeKernel(clCommandQueue(), filter.convHorizontalKernel, 2, NULL, globalWorkSizeH, filterLocalSize, 0, NULL, NULL);
	clFinish(clCommandQueue());

	clEnqueueNDRangeKernel(clCommandQueue(), filter.convVerticalKernel, 2, NULL, globalWorkSizeV, filterLocalSize, 0, NULL, NULL);
	clFinish(clCommandQueue());

	releaseGLBuffer(inImage);
	releaseGLBuffer(filter.normalsDepthBuffer);
	releaseGLBuffer(filter.tempBuffer);
	releaseGLBuffer(filter.resultBuffer);

	return filter.resultTexture;
}