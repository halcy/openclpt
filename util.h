/**
 * Some utilities useful for OpenCL hosts.
 *
 * (c) L. Diener, 2013
 */

#ifndef __UTIL_H__
#define __UTIL_H__

#include <CL/opencl.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

char *read_file(const char *filename, size_t *length);

cl_platform_id ocl_find_platform(cl_device_type device_type);
cl_device_id ocl_get_device(cl_platform_id platform, cl_device_type device_type);
cl_program ocl_build_program(cl_context context, cl_device_id device, const char* source, size_t size);
void* ocl_fetch_data(cl_command_queue command_queue, cl_mem buffer, size_t size);

const char* ocl_status_to_string(cl_int status);

float luminance(float r, float g, float b);
float write_bmp(char* output_name, int width, int height, float* r_done, float* g_done, float* b_done, float exposure);

#endif
