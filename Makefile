all:
	gcc -std=c99 -I. -L. main.c util.c -lOpenCL -o openclpt