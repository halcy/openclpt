/**
 * A simple OpenCL path tracer rendering an animation.
 * Some configuration below, for scene setup look at scene.h
 *
 * If you change things here, you might have to also change them in
 * path_tracer.cl.
 *
 * (c) L. Diener, 2013
 */

// Which OpenCL device type to use - CL_DEVICE_TYPE_GPU or CL_DEVICE_TYPE_CPU
#define OPENCL_DEVICE CL_DEVICE_TYPE_GPU

// Prefix for output file names
#define OUTPUT_PREFIX "cornell"

// How many and which (counting from 0) frames to render
#define FRAME_COUNT 90
#define FRAME_FIRST 0
#define FRAME_LAST 89
#define FRAME_STEP 1

// Image size - must match the kernel, obviously
#define IMAGE_WIDTH 400
#define IMAGE_HEIGHT 400

// Brightness scale factor to make images bright or moody
#define OUTPUT_BRIGHTNESS 7.5f
 
// Standard includes
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// OpenCL includes
#include <CL/opencl.h>

// Application includes
#include "util.h"
#include "scene.h"

// Path-trace an image and write it to disk
int main(int argc, char** argv) {
	int status;
    fprintf(stderr, "Initializing OpenCL\n");
    
	// Get a platform that hopefully has the device we want
	cl_platform_id platform = ocl_find_platform(OPENCL_DEVICE);
	
	// Get a device from the platform
	cl_device_id device = ocl_get_device(platform, OPENCL_DEVICE);
   
    // Finally, create a context and command queue for that device
	cl_context context = clCreateContext(0, 1, &device, NULL, NULL, NULL);
    cl_command_queue command_queue = clCreateCommandQueue(context, device, 0, NULL);
	
    // Load and build the OpenCL program
	fprintf(stderr, "Setting up program and kernel\n");
	size_t size;
    char* opencl_source = read_file("path_tracer.cl", &size);
    cl_program program = ocl_build_program(context, device, opencl_source, size);

    // Create kernel
    cl_kernel kernel = clCreateKernel(program, "path_tracer", &status);
	fprintf(stderr, "Kernel creation: %s\n", ocl_status_to_string(status));
	
	// Iterate over frames
	float initial_luminance = 0.0f;
	float running_luminance = 0.0f;
	for(int frame = FRAME_FIRST; frame <= FRAME_LAST; frame+=FRAME_STEP) {
		fprintf(stderr, "Rendering frame %d\n", frame);
	
		// Prepare Scene
		fprintf(stderr, "Setting up data buffers\n");
		object scene[SCENE_OBJECT_COUNT]; // SCENE_OBJECT_COUNT is defined in scene.h
		float animation_frame = ((2.0f*M_PI)/(float)FRAME_COUNT) * frame;
		prepare_scene(scene, animation_frame);
		
		// Output buffers
		int image_size = IMAGE_WIDTH * IMAGE_HEIGHT;
		cl_mem buffer_r = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * image_size, NULL, NULL);
		cl_mem buffer_g = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * image_size, NULL, NULL);
		cl_mem buffer_b = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * image_size, NULL, NULL);
		
		// Input buffer
		cl_mem buffer_scene = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(object) * SCENE_OBJECT_COUNT, scene, NULL);
		
		// Set as buffers up as arguments
		fprintf(stderr, "Setting up kernel parameters\n");
		clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &buffer_r);
		clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &buffer_g);
		clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &buffer_b);
		clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &buffer_scene);
		
		// Block until everything is ready
		clFinish(command_queue);
			
		// Execute
		fprintf(stderr, "Running kernel\n");
		status = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, (size_t*)&image_size, NULL, 0, NULL, NULL);
		fprintf(stderr, "Queueing status: %s\n", ocl_status_to_string(status));
		
		// Block until done
		clFinish(command_queue);
			
		// Fetch data
		fprintf(stderr, "Fetching results\n");
		float* r_done = (float*)ocl_fetch_data(command_queue, buffer_r, IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(float));
		float* g_done = (float*)ocl_fetch_data(command_queue, buffer_g, IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(float));
		float* b_done = (float*)ocl_fetch_data(command_queue, buffer_b, IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(float));

		// Get maximum luminance
		float max_luminance = 0.0f;
		for(int i = 0; i < IMAGE_WIDTH * IMAGE_HEIGHT; i++) {
			max_luminance = fmax(max_luminance, luminance(r_done[i], g_done[i], b_done[i]));      
		}
		
		// Set initial and running luminance
		if(initial_luminance == 0.0f) {
			initial_luminance = max_luminance;
			running_luminance = max_luminance;
			fprintf(stderr, "Initial luminance*samples: %f\n", max_luminance );
		}
		else {
			running_luminance  = (running_luminance + 1.05 * max_luminance) / 2.05f;
			fprintf(stderr, "Frame luminance*samples: %f\n", max_luminance );
			fprintf(stderr, "Running luminance*samples: %f\n", running_luminance );
		}
		
		// Tone-map and write data, first with per-frame tone mapping
		char output_name[128];
		sprintf(output_name, "perframe_%s_%05d.bmp", OUTPUT_PREFIX, frame);
		write_bmp(output_name, IMAGE_WIDTH, IMAGE_HEIGHT, r_done, g_done, b_done, max_luminance / OUTPUT_BRIGHTNESS );
		
		// Write with tone-mapping using initial frame data
		sprintf(output_name, "initial_%s_%05d.bmp", OUTPUT_PREFIX, frame);
		write_bmp(output_name, IMAGE_WIDTH, IMAGE_HEIGHT, r_done, g_done, b_done, initial_luminance / OUTPUT_BRIGHTNESS);

		// Write with tone-mapping using running exponential luminance
		sprintf(output_name, "running_%s_%05d.bmp", OUTPUT_PREFIX, frame);
		write_bmp(output_name, IMAGE_WIDTH, IMAGE_HEIGHT, r_done, g_done, b_done, running_luminance / OUTPUT_BRIGHTNESS);
		
		// If on last frame, wind down running luminance and write those frames, too
		if(frame == FRAME_LAST) {
			int bonus_frame = 0;
			while(fabs(running_luminance - initial_luminance) > OUTPUT_BRIGHTNESS * 255.0f) {
				bonus_frame++;
				fprintf(stderr, "Writing running luminance frame %d\n", frame+bonus_frame);
				running_luminance  = (running_luminance + 1.15 * max_luminance) / 2.15f;
				sprintf(output_name, "running_%s_%05d.bmp", OUTPUT_PREFIX, frame+bonus_frame);
				write_bmp(output_name, IMAGE_WIDTH, IMAGE_HEIGHT, r_done, g_done, b_done, running_luminance / OUTPUT_BRIGHTNESS);
			}
		}
		
		// Release fetch buffers and OpenCL buffers
		fprintf(stderr, "Releasing buffers\n");
		free(r_done);
		free(g_done);
		free(b_done);
		clReleaseMemObject(buffer_r);
		clReleaseMemObject(buffer_g);
		clReleaseMemObject(buffer_b);
		clReleaseMemObject(buffer_scene);
	}
	
	// Release OpenCL ressources
	fprintf(stderr, "Shutting down\n");
	clReleaseKernel(kernel); 
    clReleaseProgram(program);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);
	
	// Done
	fprintf(stderr, "Done\n");
    exit(0);
}
