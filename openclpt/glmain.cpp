/**
 * A simple, but neat path tracer.
 *
 * L. Diener, 2014
 */

// Window settings
#define WINDOW_TITLE "OpenCLPT"

int WINDOW_WIDTH = 512;
int WINDOW_HEIGHT = 512;

// #define FULLSCREEN

// OpenGL Debug Mode toggle
#define DEBUG

// "THIS FUNCTION IS UNSAFE"
#pragma warning(disable: 4996)

// Include files with neat things.
#include "glhelpers.h"
#include "obj_loader.h"
#include "speedup_grid.h"
#include "filter.h"
#include "bmp_handler.h"

// Vector tools
#include "Vector.h"

#include <math.h>
#include <time.h>

// Include CL/GL interop tools.
#include "clgl.h"

// Data structures that can be sent over to OpenCL
#include "data_structures.h"

// Data like shaders and such
#include "data.h"

// Vertex array object
GLuint vertexArray;

// User brightness scale
float userScale = 1.0f;

// Switch to run or not run parts
bool runTracer = true;
bool runFilter = false;
bool writeBitmap = false;
bool cameraLock = true;
bool cameraChanged = false;
bool recalcFilter = true;

// For buffer clears
float* zeroPixels;

// Forward declaration
void draw();
void update();
void updateI(int);
void handleKeypress(unsigned char k, int x, int y);

// Per pixel RNG texture
cl_mem randBuffer;

// How many samples per kernel invocation
int samplesPerCall = 1;

//////////////////////////////// INIT FUNCTIONS ////////////////////////////////

// Make a window for doing OpenGL
void makeWindow(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitContextVersion(4, 2);
	#ifdef DEBUG
		glutInitContextFlags(GLUT_CORE_PROFILE | GLUT_DEBUG);
	#else
		glutInitContextFlags(GLUT_COMPATIBILITY_PROFILE);
	#endif

	// Open window / full screen
	#ifdef FULLSCREEN
		char modeString[255] = "";
		sprintf(modeString, "%dx%d:24", WINDOW_WIDTH, WINDOW_HEIGHT);
		glutGameModeString(modeString);
		glutEnterGameMode();
	#else
		glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
		glutCreateWindow(WINDOW_TITLE);
	#endif

	// Use GLEW, and make sure it imports EVERYTHING
	glewExperimental = GL_TRUE;
	glewInit();

	#ifdef DEBUG
		registerGlDebugLogger(GL_DEBUG_SEVERITY_MEDIUM);
	#endif

	// Set up OpenGL features
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0f);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
}

// Functions that glut calls for us
void setupGlutCallbacks() {
	glutDisplayFunc(draw);
	glutTimerFunc(15, updateI, 0);
	glutKeyboardFunc(handleKeypress);
}

// Initialize buffers, textures, etc.
void initObjects() {
	// Create a VAO and bind it
	glGenVertexArrays(1, &vertexArray);
	glBindVertexArray(vertexArray);

	// Prepare a screen quad to render postprocessed things.
	vec3_t quadData[] = {
		{-1.0f, -1.0f, 0.0f, 1.0f},
		{ 1.0f, -1.0f, 0.0f, 1.0f},
		{ 1.0f,  1.0f, 0.0f, 1.0f},
		{-1.0f,  1.0f, 0.0f, 1.0f}
	};
	GLuint quadElements[] = {0, 1, 3, 1, 2, 3};

	screenQuad.vertexBuffer = makeBO(
		GL_ARRAY_BUFFER,
		quadData,
		sizeof(vec3_t) * 4,
		GL_STATIC_DRAW
	);
	screenQuad.elementBuffer = makeBO(
		GL_ELEMENT_ARRAY_BUFFER,
		quadElements,
		sizeof(GLuint) * 6,
		GL_STATIC_DRAW
	);

	// Create textures and share them with OpenCL
	zeroPixels = (float*)malloc(WINDOW_WIDTH * WINDOW_HEIGHT * 4 * sizeof(float));
	for(int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT * 4; i++) {
		zeroPixels[i] = 0.0f;
	}
	for(int i = 0; i < 2; i++) {
		displayShader.texture[i] = makeTextureBuffer(WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGBA, GL_RGBA32F);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RGBA, GL_FLOAT, zeroPixels);
		tracer.renderBuffer[i] = clCreateFromGLTexture2D(clContext(), CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, displayShader.texture[i], 0);
	}

	// Start writing to buffer 0
	tracer.targetBufferIndex = 0;
	tracer.sampleCount = 0.0f;
}

// Create an OpenCL program and kernels
void initPrograms() {
	char defines[1024];
	sprintf(defines, "-D IMAGE_WIDTH=%d -D IMAGE_HEIGHT=%d -D GRID_RADIUS=%d -D GRID_RES=%f -D SAMPLES=%d", WINDOW_WIDTH, WINDOW_HEIGHT, scene.gridRadius, scene.gridRes, samplesPerCall);
	tracer.program = clProgramFromFile("tracer.cl", defines);
	tracer.renderKernel = clCreateKernel(tracer.program, "PathTracer", NULL);
	tracer.clearKernel = clCreateKernel(tracer.program, "ClearImage", NULL);
	prepareFilter(tracer.program, WINDOW_WIDTH, WINDOW_HEIGHT);
}

// Initialize shaders.
void initShaders() {
	// Load a simple display-that-thing Vertex/Fragment shader
	GLuint vertexShader = loadShader(GL_VERTEX_SHADER, "quad.vert");
	GLuint fragmentShader = loadShader(GL_FRAGMENT_SHADER, "simple.frag");
	displayShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);

	// Get locations of attributes and uniforms used inside.
	displayShader.vertexPosition = glGetAttribLocation(displayShader.shaderProgram, "vertex");
	displayShader.textureLocation = glGetUniformLocation(displayShader.shaderProgram, "displayTexture");
	displayShader.iterCount =  glGetUniformLocation(displayShader.shaderProgram, "iterCount");
	displayShader.userScale =  glGetUniformLocation(displayShader.shaderProgram, "userScale");
	
	// Bind output variables
	glBindFragDataLocation(displayShader.shaderProgram, 0, "outColor");
}

//////////////////////// UPDATING AND DRAWING ////////////////////////////////////////

// Update status variables.
// Called every 15ms, unless the PC is too slow.
void updateI(int) { update(); }

void update() {
	// Redraw screen now, please, and call again in 15ms.
	glutPostRedisplay();
	glutTimerFunc(15, updateI, 0);

	if(!cameraLock) {
		HWND window = GetActiveWindow();
		POINT p;
		GetCursorPos(&p);
		ScreenToClient(window, &p);
		glutWarpPointer(WINDOW_WIDTH / 2, WINDOW_HEIGHT / 2);

		float angleX = (p.x - (WINDOW_WIDTH / 2)) * CAM_ROTSPEED;
		Quaternion rotX = RotationQuaternion(-angleX, MakeVector(0, 1, 0));
		camera.front = TransformVector(RotationMatrixFromQuaternion(rotX), camera.front);

		float angleY = (p.y - (WINDOW_HEIGHT / 2)) * CAM_ROTSPEED;
		camera.elevation += angleY;
		camera.elevation = max(-0.9f, min(camera.elevation, 0.9f));

		if(fabs(angleX) + fabs(angleY) > 0.001f) {
			cameraChanged = true;
		}

		if(GetAsyncKeyState('W') != 0) {
			camera.pos = VectorAdd(camera.pos, VectorMul(camera.front, CAM_MOVESPEED));
			cameraChanged = true;
		}

		if(GetAsyncKeyState('S') != 0) {
			camera.pos = VectorAdd(camera.pos, VectorMul(camera.front, -CAM_MOVESPEED));
			cameraChanged = true;
		}

		if(GetAsyncKeyState('E') != 0) {
			camera.pos = VectorAdd(camera.pos, VectorMul(camera.up, CAM_MOVESPEED));
			cameraChanged = true;
		}

		if(GetAsyncKeyState('Q') != 0) {
			camera.pos = VectorAdd(camera.pos, VectorMul(camera.up, -CAM_MOVESPEED));
			cameraChanged = true;
		}

		if(GetAsyncKeyState('D') != 0) {
			camera.pos = VectorAdd(camera.pos, VectorMul(VectorCross(camera.front, camera.up), CAM_MOVESPEED));
			cameraChanged = true;
		}

		if(GetAsyncKeyState('A') != 0) {
			camera.pos = VectorAdd(camera.pos, VectorMul(VectorCross(camera.front, camera.up), -CAM_MOVESPEED));
			cameraChanged = true;
		}
	}
}

// Scene loader
void loadScene(char* sceneFile, char* matFile, int gridExtent) {
	// Load actual scene and materials
	scene.scene = loadObj(sceneFile, &scene.triCount, matFile);

	// Set up grid
	scene.gridRadius = gridExtent;
	vec3_t maxVal = vec3(-1000000.0f, -1000000.0f, -1000000.0f);
	vec3_t minVal = vec3( 1000000.0f,  1000000.0f,  1000000.0f);
	for(int i = 0; i < scene.triCount; i++) {
		maxVal = vec3(
			max(maxVal.x, scene.scene[i].object_data.p1.x),
			max(maxVal.y, scene.scene[i].object_data.p1.y),
			max(maxVal.z, scene.scene[i].object_data.p1.z)
		);
		maxVal = vec3(
			max(maxVal.x, scene.scene[i].object_data.p2.x),
			max(maxVal.y, scene.scene[i].object_data.p2.y),
			max(maxVal.z, scene.scene[i].object_data.p2.z)
		);
		maxVal = vec3(
			max(maxVal.x, scene.scene[i].object_data.p2.x),
			max(maxVal.y, scene.scene[i].object_data.p2.y),
			max(maxVal.z, scene.scene[i].object_data.p2.z)
		);
		minVal = vec3(
			min(minVal.x, scene.scene[i].object_data.p1.x),
			min(minVal.y, scene.scene[i].object_data.p1.y),
			min(minVal.z, scene.scene[i].object_data.p1.z)
		);
		minVal = vec3(
			min(minVal.x, scene.scene[i].object_data.p2.x),
			min(minVal.y, scene.scene[i].object_data.p2.y),
			min(minVal.z, scene.scene[i].object_data.p2.z)
		);
		minVal = vec3(
			min(minVal.x, scene.scene[i].object_data.p2.x),
			min(minVal.y, scene.scene[i].object_data.p2.y),
			min(minVal.z, scene.scene[i].object_data.p2.z)
		);
	}
	float maxTotal = 0.0f;
	maxTotal = max(maxTotal, fabs(maxVal.x));
	maxTotal = max(maxTotal, fabs(maxVal.y));
	maxTotal = max(maxTotal, fabs(maxVal.z));
	maxTotal = max(maxTotal, fabs(minVal.x));
	maxTotal = max(maxTotal, fabs(minVal.y));
	maxTotal = max(maxTotal, fabs(minVal.z));
	scene.gridRes = (maxTotal / (((float)scene.gridRadius))) * 1.03f;

	printf("Calculating grid with extent %d, resolution %f\n", scene.gridRadius, scene.gridRes);
	scene.gridSceneSize = createGrid(scene.scene, scene.triCount, scene.gridRadius, scene.gridRes, false, &scene.gridScene, &scene.gridOffsets);
	printf("Grid done, %d entries\n", scene.gridSceneSize);

	scene.onDevice = false;
}

// Draw the scene to the screen
void draw() {
	//////////////////////// PART 1: RENDERIING ////////////////////////////////

	if(runTracer) {
		// Grab buffers for OpenCL
		acquireGLBuffer(tracer.renderBuffer[0]);
		acquireGLBuffer(tracer.renderBuffer[1]);

		// Prepare to run some kernels
		cl_int numPixels = WINDOW_WIDTH * WINDOW_HEIGHT;

		cl_uint workSize[3] = {numPixels, 0, 0};
		cl_uint workgroupSize[3] = {256, 0, 0};

		// Send scene
		if(scene.onDevice == false) {
			scene.sceneBuffer = clCreateBuffer(
				clContext(), 
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
				sizeof(object) * scene.gridSceneSize, 
				scene.gridScene, 
				NULL
			);
			scene.offsetBuffer = clCreateBuffer(
				clContext(), 
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
				sizeof(int) * (2 * scene.gridRadius) * (2 * scene.gridRadius) * (2 * scene.gridRadius), 
				scene.gridOffsets, 
				NULL
			);
			scene.onDevice = true;

			// Init rng
			int* randomData = (int*)malloc(sizeof(int) * WINDOW_WIDTH * WINDOW_HEIGHT);
			for(int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; i++) {
				randomData[i] = rand();
			}
			randBuffer = clCreateBuffer(
				clContext(), 
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
				sizeof(int) * WINDOW_WIDTH * WINDOW_HEIGHT, 
				randomData, 
				NULL
			);
			free(randomData);

			// Set up static arguments
			clSetKernelArg(tracer.renderKernel, 2, sizeof(cl_mem), &randBuffer);
			clSetKernelArg(tracer.renderKernel, 3, sizeof(cl_mem), &scene.sceneBuffer);
			clSetKernelArg(tracer.renderKernel, 4, sizeof(cl_mem), &scene.offsetBuffer);
		}

		// Potentially clear input buffer
		if(cameraChanged) {
			clSetKernelArg(tracer.clearKernel, 0, sizeof(cl_mem), &tracer.renderBuffer[1 - tracer.targetBufferIndex]);
			clRunKernel(tracer.clearKernel, workSize, workgroupSize);
			cameraChanged = false;
			tracer.sampleCount = 0;
			recalcFilter = true;
		}

		// Potentially refill RNG
		if((int)tracer.sampleCount % 10 == 0) {
			int* randomData = (int*)malloc(sizeof(int) * WINDOW_WIDTH * WINDOW_HEIGHT);
			for(int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; i++) {
				randomData[i] = rand();
			}
			clEnqueueWriteBuffer(clCommandQueue(), randBuffer, true, 0, sizeof(int) * WINDOW_WIDTH * WINDOW_HEIGHT, randomData, 0, 0, 0);
			clFinish(clCommandQueue());
			free(randomData);
		}

		// Set up arguments
		clSetKernelArg(tracer.renderKernel, 0, sizeof(cl_mem), &tracer.renderBuffer[1 - tracer.targetBufferIndex]);
		clSetKernelArg(tracer.renderKernel, 1, sizeof(cl_mem), &tracer.renderBuffer[tracer.targetBufferIndex]);

		float cameraPos[4];
		cameraPos[0] = -camera.pos.x;
		cameraPos[1] = -camera.pos.y;
		cameraPos[2] = -camera.pos.z;
		cameraPos[3] = 0.0f;
		clSetKernelArg(tracer.renderKernel, 5, sizeof(cl_float4), &cameraPos);

		Vector vectorUp = MakeVector(0, 1, 0);
		Quaternion rotZ = RotationQuaternion(-camera.elevation, VectorNorm(VectorCross(VectorNorm(camera.front), vectorUp)));

		float cameraUp[4];
		Vector cameraUpV = VectorNorm(TransformVector(RotationMatrixFromQuaternion(rotZ), vectorUp));
		cameraUp[0] = cameraUpV.x;
		cameraUp[1] = cameraUpV.y;
		cameraUp[2] = cameraUpV.z;
		cameraUp[3] = 0.0f;
		clSetKernelArg(tracer.renderKernel, 6, sizeof(cl_float4), &cameraUp);

		float cameraFront[4];
		Vector cameraFrontV = VectorNorm(TransformVector(RotationMatrixFromQuaternion(rotZ), camera.front));
		cameraFront[0] = cameraFrontV.x;
		cameraFront[1] = cameraFrontV.y;
		cameraFront[2] = cameraFrontV.z;
		cameraFront[3] = 0.0f;
		clSetKernelArg(tracer.renderKernel, 7, sizeof(cl_float4), &cameraFront);

		// Trace an iteration
		clRunKernel(tracer.renderKernel, workSize, workgroupSize);

		// Sample has been registered
		tracer.sampleCount = tracer.sampleCount + (float)samplesPerCall;

		// Release buffers back to OpenGL
		releaseGLBuffer(tracer.renderBuffer[0]);
		releaseGLBuffer(tracer.renderBuffer[1]);
	}

	//////////////////////// PART 1.5: FILTER ///////////////////////////////

	GLuint displayTexture;
	if(runFilter) {
		if(recalcFilter) {
			float cameraPos[4];
			cameraPos[0] = -camera.pos.x;
			cameraPos[1] = -camera.pos.y;
			cameraPos[2] = -camera.pos.z;
			cameraPos[3] = 0.0f;

			Vector vectorUp = MakeVector(0, 1, 0);
			Quaternion rotZ = RotationQuaternion(-camera.elevation, VectorNorm(VectorCross(VectorNorm(camera.front), vectorUp)));

			float cameraUp[4];
			Vector cameraUpV = VectorNorm(TransformVector(RotationMatrixFromQuaternion(rotZ), vectorUp));
			cameraUp[0] = cameraUpV.x;
			cameraUp[1] = cameraUpV.y;
			cameraUp[2] = cameraUpV.z;
			cameraUp[3] = 0.0f;

			float cameraFront[4];
			Vector cameraFrontV = VectorNorm(TransformVector(RotationMatrixFromQuaternion(rotZ), camera.front));
			cameraFront[0] = cameraFrontV.x;
			cameraFront[1] = cameraFrontV.y;
			cameraFront[2] = cameraFrontV.z;
			cameraFront[3] = 0.0f;

			// Prepare filter normaldepth
			displayTexture = prepareNormalDepth(randBuffer, scene.sceneBuffer, scene.offsetBuffer, cameraPos, cameraUp, cameraFront);
			recalcFilter = false;
		}

		displayTexture = runFilterKernels(tracer.renderBuffer[tracer.targetBufferIndex]);
	}
	else {
		displayTexture = displayShader.texture[tracer.targetBufferIndex];
	}

	//////////////////////// PART 2: DISPLAY ////////////////////////////////

	// Clear everything first thing.
	glClearColor(0.0f, 1.0f, 0.0f, 1.0f);

	// Compose shader
	glUseProgram(displayShader.shaderProgram);

	// Bind and set textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, displayTexture);
	glUniform1i(displayShader.textureLocation, 0);
	
	// Send uniforms
	glUniform1f(displayShader.iterCount, tracer.sampleCount);
	glUniform1f(displayShader.userScale, userScale);
	
	// Draw a quad
    glDisable(GL_DEPTH_TEST);
	glBindBuffer(GL_ARRAY_BUFFER, screenQuad.vertexBuffer);
	glVertexAttribPointer(
		displayShader.vertexPosition,
		3,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 4,
		(void*)0
	);
	glEnableVertexAttribArray(displayShader.vertexPosition);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, screenQuad.elementBuffer);
	glDrawElements(
		GL_TRIANGLES,
		6,
		GL_UNSIGNED_INT,
		(void*)0
	);
	glEnable(GL_DEPTH_TEST);

	//////////////////////// PART 3: WRITE ////////////////////////////////
	if(writeBitmap == true) {
		writeBitmap = false;
		char fileName[1024];
		sprintf(fileName, "tracer_out_%d_%d.bmp", (int)tracer.sampleCount, time(NULL));
		bmp_init(fileName, WINDOW_WIDTH, WINDOW_HEIGHT);
		glBindTexture(GL_TEXTURE_2D, displayTexture);
		float* pixelData = (float*)malloc(sizeof(float) * 4 * WINDOW_WIDTH * WINDOW_HEIGHT);
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, pixelData);
		for(int i = 0; i < 4 * WINDOW_WIDTH * WINDOW_HEIGHT; i += 4) {
			int r = (int)max(0.0f, min(((pixelData[i] / userScale) / tracer.sampleCount) * 255.0f, 255.0f));
			int g = (int)max(0.0f, min(((pixelData[i + 1] / userScale) / tracer.sampleCount) * 255.0f, 255.0f));
			int b = (int)max(0.0f, min(((pixelData[i + 2] / userScale) / tracer.sampleCount) * 255.0f, 255.0f));
			bmp_pixel(r, g, b);
		}
		bmp_close();
		printf("Wrote current buffer to %s\n", fileName);
	}

	if(runTracer == true) {
		// Switch OpenCL target buffer
		tracer.targetBufferIndex = 1 - tracer.targetBufferIndex;
	}

	// Switch drawing area and displayed area.
	glutSwapBuffers();
}

// Key press handler
// Mostly here to allow us to break on escape
void handleKeypress(unsigned char k, int x, int y) {
	FILE* camFile;
	switch(k) {
		case 27: // Escape -> die.
			exit(0);
		break;

		case 'f':
			userScale += 0.001f;
		break;

		case 'g':
			userScale -= 0.001f;
		break;

		case '1':
			userScale += 0.01f;
		break;

		case '2':
			userScale -= 0.01f;
		break;

		case '3':
			userScale += 0.1f;
		break;

		case '4':
			userScale -= 0.1f;
		break;

		case '5':
			userScale += 1.0f;
		break;

		case '6':
			userScale -= 1.0f;
		break;

		case '7':
			userScale += 10.0f;
		break;

		case '8':
			userScale -= 10.0f;
		break;

		case '9':
			userScale = 1.0f;
		break;

		case 'y':
			cameraLock = true;
			runTracer = false;
		break;

		case 'x':
			runTracer = true;
		break;

		case 'c':
			runFilter = false;
		break;

		case 'v':
			cameraLock = true;
			runFilter = true;
		break;

		case 'b':
			writeBitmap = true;
		break;

		case 'l':
			runTracer = true;
			cameraLock = !cameraLock;
			if(cameraLock) {
				glutSetCursor(GLUT_CURSOR_INHERIT);
			}
			else {
				glutSetCursor(GLUT_CURSOR_NONE);
			}
		break;

		case 'k':
			camFile = fopen("camera.cam", "w");
			fprintf(camFile, "%f %f %f %f %f %f %f", camera.pos.x, camera.pos.y, camera.pos.z, camera.front.x, camera.front.y, camera.front.z, camera.elevation);
			fclose(camFile);
			printf("Camera data dumped to camera.cam\n");
		break;

		case 'p': // Neat for debugging. Wait a second on 'p'.
			Sleep(1000);
		break;
	}
}

// Set things up and run
int main(int argc, char** argv) {
	WINDOW_WIDTH = atoi(argv[1]);
	WINDOW_HEIGHT = atoi(argv[2]);
	int gridExtent = atoi(argv[3]);
	char* sceneFile = argv[4];
	char* matFile = argv[5];
	samplesPerCall = atoi(argv[6]);
	char* camFileName = argv[7];

	if(argc != 8) {
		printf("Missing arguments. See README for details.\n");
	}

	printf(
		"Initializing:\n%dx%d window\nGrid extent %d\nScene %s\nMaterial library %s\nSamples per call %d\nCamera %s\n", 
		WINDOW_WIDTH, 
		WINDOW_HEIGHT,
		gridExtent,
		sceneFile, 
		matFile, 
		samplesPerCall, 
		camFileName
	);

	camera.up = MakeVector(0, 1, 0);
	FILE* camFile = fopen(camFileName, "r");
	fscanf(camFile, "%f %f %f %f %f %f %f", &camera.pos.x, &camera.pos.y, &camera.pos.z, &camera.front.x, &camera.front.y, &camera.front.z, &camera.elevation);
	fclose(camFile);
	printf("Initial camera pos: %f %f %f\n", camera.pos.x, camera.pos.y, camera.pos.z);
	printf("Initial camera front: %f %f %f\n", camera.front.x, camera.front.y, camera.front.z);
	printf("Initial camera elevation: %f\n", camera.elevation);

	loadScene(sceneFile, matFile, gridExtent);
	makeWindow(argc, argv);
	acquireSharedOpenCLContext();
	initObjects();
	initPrograms();
	initShaders();
	setupGlutCallbacks();
	glutMainLoop();
	releaseSharedOpenCLContext();
}