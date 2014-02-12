OpenCLPT
==========

A simple OpenCL path tracer with interactive preview.

Command line options (Mandatory, see also included .bat files):

openclpt.exe 
	<window width> \
	<window height> \
	<speedup grid extent> \
	<scene file> \
	<material file> \
	<samples per kernel invocation> \
	<camera initialization file>

Controls: 
* Tracing:
	* Y to pause path tracing
	* X to resume path tracing
* Movement:
	* L to lock/unlock camera
	* Mouse to look
	* WASD to move around
	* Q/E to go up/down
* Preview and Output:
	* f/g, 1/2 through 7/8: Decrease / Increase brightness slightly -> strongly
	* B to dump current frame buffer to bitmap file
	* K to dump current camera position and orientation to file
* Post-filtering:
	* V: Activate bilateral filter
	* C: Deactivate bilateral filter
* Misc:
	* P to pause everything for a second
	* Escape to quit

To make, open the .sln file in Visual Studio and compile, everything 
that is required should be included in the repository. Currently,
only windows is supported, but most of the code (everything that is
not user input handling) is platform independent and should work on
any platform with a new enough OpenGL.

![Screenshot](https://github.com/halcy/openclpt/blob/master/coolbox.png?raw=true)

The path tracer is relatively simple and straightforward, paths are traced
using standard monte-carlo path tracing with russian roulette termination.
Tracing is accelerated using a regular grid data structure pre-computed on 
the CPU (Precomputation should only take a few seconds even for complex
scenes).

The input file format for scenes and materials is a variant of the obj/mtl
format: For the scene itself, the format is a subset of obj - only triangles
with vertex normals are supported. Material names must be numeric. The material
format is as follows:

m <material number>
a <albedo r> <albedo g> <albedo b>
e <emissive r> <emissive g> <emissive b>
p <specularity (0 to 1)> <transparency (0 to 1)> <specular quality>

The albedo is simply the colour of the material, the emissive colour is the 
colour of light emitted. Specularity gives how much of the material is 
specular / transmissive. Transparency gives how much of the specular
part is transmissive. The specular quality controls the glossiness of
the material - higher means more perfect reflections. See the included
material library for some examples.
