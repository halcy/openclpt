openclpt
========

A simple OpenCL path tracer with spheres and mixed glossy/specular + diffuse + emission material.

If run as set up, it will render a small 90 frame animation, using the GPU, at 400x400 with
30kspp. On my GPU (Reasonably new NVidia GeForce), it takes about 8 to 12 minutes to render
each frame. Note that during the time it renders, it pretty much locks down your computer 
entirely, since your GPU cannot do things like "update what is on screen" while it runs. If
you just want to try it out, I recommend setting the samples per pixel to maybe 300, and 
rendering a single frame. You can also try using the CPU, to see just how much slower it is
(It is orders of magnitude slower). While the code probably isn't very optimal or even very
nice, it is commented and simple, unlike most OpenCL "example" monstrosities - you might be
able to just read through this and kind of understand what is going on!

Here is what the output is like, if you don't want to wait: 
 
 http://aka-san.halcy.de/tracer/tracer.mp4

To make:

> $ make

You will need OpenCL headers, obviously, but that is it! If, you, like me, only have Windows
on the machine that has the powerful GPU, that's not a problem, this code compiles and runs
fine under MinGW / msys! Just copy the libOpenCL.a and headers from the windows_support folder 
into the root folder and make sure you have OpenCL installed, it should(tm) work fine! (It
works fine for me, at least - I have no machine with AMD GPU to test on, though.)

To run:

> $ ./clpt

p.s. for fun, try playing around with the PRNG in path_tracer.cl, the code is very sensitive
to strangely behaving PRNGs! (If you come up with anything fun, let me know!)
