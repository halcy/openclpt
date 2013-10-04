/**
 * Some utilities useful for OpenCL hosts or path tracers
 *
 * (c) L. Diener, 2013
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"
#include "bmp_handler.h"

// Read a file
char *read_file(const char* name, size_t* length) {
    FILE *f = fopen(name, "r");

	// Get length
    fseek(f, 0, SEEK_END);
	*length = ftell(f);
    fseek(f, 0, SEEK_SET);

	// Read
    char *buffer = (char*)malloc(*length+1);
	*length = fread(buffer, 1, *length, f);
	fclose(f);
	
	// Be paranoid
    buffer[*length] = '\0';

    return buffer;
}

// Get an OpenCL platform id
cl_platform_id ocl_find_platform(cl_device_type device_type) {
    char name[128];

    // Get OpenCL platform count
	cl_uint platform_count;
    clGetPlatformIDs(0, NULL, &platform_count);
	if(platform_count == 0) {
		fprintf(stderr, "No OpenCL platforms - aborting!");
        exit(-1);
    }
	
    // Grab platform info
	cl_platform_id platform_ids[platform_count];
	clGetPlatformIDs(platform_count, platform_ids, NULL);
	
	// Assumption: CPU compute device is device 0.
	// Probably not a good assumption to make. Probably fine anyways.
	if(device_type == CL_DEVICE_TYPE_CPU) {
		clGetPlatformInfo (platform_ids[0], CL_PLATFORM_NAME, 128, &name, NULL);
		fprintf(stderr, "Using platform %d: %s\n", 0, name);
		return platform_ids[0];
	}
	
	// Find the first device after device 0 that has "NVIDIA" or "ATI" 
	// or "AMD" in the name. Failing that, return the last device.
	for(int i = 1; i < platform_count; i++) {
		clGetPlatformInfo (platform_ids[i], CL_PLATFORM_NAME, 128, &name, NULL);
        if(
			strstr(name, "NVIDIA") != NULL || 
			strstr(name, "ATI") != NULL || 
			strstr(name, "AMD") != NULL
		) {
			fprintf(stderr, "Using platform %d: %s\n", i, name);
			return platform_ids[i];
		}
    }
	fprintf(stderr, "Using platform %d: %s\n", platform_count-1, name);
	return platform_ids[platform_count-1];
}

// Given a platform and device type, get a device and return a device context.
cl_device_id ocl_get_device(cl_platform_id platform, cl_device_type device_type) {
	// Get the number of devices
	cl_uint device_count = 0;
	clGetDeviceIDs(platform, device_type, 0, NULL, &device_count);
    
	// Abort if none
	if(device_count == 0) {
		fprintf(stderr, "No fitting OpenCL device found - aborting!");
        exit(-1);
    }

    // Get device list
    cl_device_id devices[device_count];
	clGetDeviceIDs(platform, device_type, device_count, devices, NULL);
 
	// Print some data about device 0
	char name[128];
	clGetDeviceInfo (devices[0], CL_DEVICE_NAME, 128, &name, NULL);
	printf("Using device 0: %s\n", name);
	
    cl_ulong mem_size;
	clGetDeviceInfo (devices[0], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
	printf(" Global memory: %u kb\n", (unsigned int)(mem_size/1024));
	
	clGetDeviceInfo (devices[0], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
	printf(" Local memory: %u kb\n", (unsigned int)(mem_size/1024));
	
	clGetDeviceInfo (devices[0], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem_size), &mem_size, NULL);
	printf(" Constant memory: %u kb\n",  (unsigned int)(mem_size/1024));
	
	size_t workgroup_sizes[3];
	clGetDeviceInfo (devices[0], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * 3, &workgroup_sizes, NULL);
	printf(" Workgroup sizes max: (%u, %u, %u)\n", workgroup_sizes[0], workgroup_sizes[1], workgroup_sizes[2]);

	// Return the id for device 0
    return devices[0];
}

// Create and compile a program and return the compiled program,
// given context, device and program source.
cl_program ocl_build_program(cl_context context, cl_device_id device, const char* source, size_t size) {
	// Create
	cl_program program = clCreateProgramWithSource(context, 1, &source, &size, NULL);
	
	// Build
    cl_int status = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    printf("Build result: %s\n", ocl_status_to_string(status));
	
	// Print build log
	size_t log_size;
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
	
	char build_log [log_size+1];
	clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL);
	
	build_log[log_size] = '\0';
	printf("Build log:\n------------------------\n%s\n------------------------\n", build_log);

	return program;
}

// Fetch some data from an OpenCL buffer back to host memory
void* ocl_fetch_data(cl_command_queue command_queue, cl_mem buffer, size_t size) {
	// Allocate and fetch
	void* data = malloc(size);
	clEnqueueReadBuffer(
		command_queue,
		buffer,
		CL_TRUE,
		0,
		size,
		data,
		0,
		NULL,
		NULL
	);
	
	// Block until done and return
	clFinish(command_queue);
	return data;
}

// Translate an OpenCL return code to a string
const char* ocl_status_to_string(cl_int status) {
    switch (status) {
        case CL_SUCCESS:                            return (const char*)"Success!"; 
        case CL_DEVICE_NOT_FOUND:                   return (const char*)"Device not found.";
        case CL_DEVICE_NOT_AVAILABLE:               return (const char*)"Device not available";
        case CL_COMPILER_NOT_AVAILABLE:             return (const char*)"Compiler not available";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return (const char*)"Memory object allocation failure";
        case CL_OUT_OF_RESOURCES:                   return (const char*)"Out of resources";
        case CL_OUT_OF_HOST_MEMORY:                 return (const char*)"Out of host memory";
        case CL_PROFILING_INFO_NOT_AVAILABLE:       return (const char*)"Profiling information not available";
        case CL_MEM_COPY_OVERLAP:                   return (const char*)"Memory copy overlap";
        case CL_IMAGE_FORMAT_MISMATCH:              return (const char*)"Image format mismatch";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return (const char*)"Image format not supported";
        case CL_BUILD_PROGRAM_FAILURE:              return (const char*)"Program build failure";
        case CL_MAP_FAILURE:                        return (const char*)"Map failure";
        case CL_INVALID_VALUE:                      return (const char*)"Invalid value";
        case CL_INVALID_DEVICE_TYPE:                return (const char*)"Invalid device type";
        case CL_INVALID_PLATFORM:                   return (const char*)"Invalid platform";
        case CL_INVALID_DEVICE:                     return (const char*)"Invalid device";
        case CL_INVALID_CONTEXT:                    return (const char*)"Invalid context";
        case CL_INVALID_QUEUE_PROPERTIES:           return (const char*)"Invalid queue properties";
        case CL_INVALID_COMMAND_QUEUE:              return (const char*)"Invalid command queue";
        case CL_INVALID_HOST_PTR:                   return (const char*)"Invalid host pointer";
        case CL_INVALID_MEM_OBJECT:                 return (const char*)"Invalid memory object";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return (const char*)"Invalid image format descriptor";
        case CL_INVALID_IMAGE_SIZE:                 return (const char*)"Invalid image size";
        case CL_INVALID_SAMPLER:                    return (const char*)"Invalid sampler";
        case CL_INVALID_BINARY:                     return (const char*)"Invalid binary";
        case CL_INVALID_BUILD_OPTIONS:              return (const char*)"Invalid build options";
        case CL_INVALID_PROGRAM:                    return (const char*)"Invalid program";
        case CL_INVALID_PROGRAM_EXECUTABLE:         return (const char*)"Invalid program executable";
        case CL_INVALID_KERNEL_NAME:                return (const char*)"Invalid kernel name";
        case CL_INVALID_KERNEL_DEFINITION:          return (const char*)"Invalid kernel definition";
        case CL_INVALID_KERNEL:                     return (const char*)"Invalid kernel";
        case CL_INVALID_ARG_INDEX:                  return (const char*)"Invalid argument index";
        case CL_INVALID_ARG_VALUE:                  return (const char*)"Invalid argument value";
        case CL_INVALID_ARG_SIZE:                   return (const char*)"Invalid argument size";
        case CL_INVALID_KERNEL_ARGS:                return (const char*)"Invalid kernel arguments";
        case CL_INVALID_WORK_DIMENSION:             return (const char*)"Invalid work dimension";
        case CL_INVALID_WORK_GROUP_SIZE:            return (const char*)"Invalid work group size";
        case CL_INVALID_WORK_ITEM_SIZE:             return (const char*)"Invalid work item size";
        case CL_INVALID_GLOBAL_OFFSET:              return (const char*)"Invalid global offset";
        case CL_INVALID_EVENT_WAIT_LIST:            return (const char*)"Invalid event wait list";
        case CL_INVALID_EVENT:                      return (const char*)"Invalid event";
        case CL_INVALID_OPERATION:                  return (const char*)"Invalid operation";
        case CL_INVALID_GL_OBJECT:                  return (const char*)"Invalid OpenGL object";
        case CL_INVALID_BUFFER_SIZE:                return (const char*)"Invalid buffer size";
        case CL_INVALID_MIP_LEVEL:                  return (const char*)"Invalid mip-map level";
        default: 									return (const char*)"Unknown";
    }
}

// Get luminance from RGB
float luminance(float r, float g, float b) {
    return 0.2126f * r + 0.7152f * g + 0.0722f * b;
}

// "Tone-map" and write an image
float write_bmp(char* output_name, int width, int height, float* r_done, float* g_done, float* b_done, float exposure) {
	bmp_init(output_name, width, height );
	for(int i = 0; i < width * height; i++) {
		float r = pow(fmax(0.0f, fmin(1.0f, r_done[i] / exposure)), 1.0f/2.2f) * 255.0f;
		float g = pow(fmax(0.0f, fmin(1.0f, g_done[i] / exposure)), 1.0f/2.2f) * 255.0f;
		float b = pow(fmax(0.0f, fmin(1.0f, b_done[i] / exposure)), 1.0f/2.2f) * 255.0f;
		bmp_pixel((int)r, (int)g, (int)b);
	}
	bmp_close();
}