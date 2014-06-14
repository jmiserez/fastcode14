fastcode14
==========

Fast Exposure Fusion implementation

- Conversion of images to uncompressed TIFF:
  $ convert -compress none file.jpg
  $ convert -resize 10% -compress none file.jpg
 
Directories
==============================
/src:
  /01ref_matlab: 1:1 matlab implementation
  /consolidated: Different versions of the source code for benchmarking
  /driver: Benchmarking support code
  /perf: PMU support code, cost model support code
  /standalone: completely standalone version of the code, with separate Makefile
  /unused: unused and unsupported code
/reference: contains the original reference Matlab code
/testdata: contains sample input and output images for verification
/utils: contains dynamorio code to measure ILP


Running benchmarks
==================
fusion_driver.py configures, builds and runs different versions 
of code for benchmarking.

cd src/
./fusion_driver.py

Options can be changed in the file fusion_driver.py, 
all Makefiles and binaries are generated automatically.


The optimized versions of the source code are in the 
src/consolidated/ folder, which are then combined using the
appropriate Make.config files.

Standalone version with debug outputs
=====================================

/src/standalone contains a standalone version of the program, it can be run
with:

cd src/standalone/
make
./run_example.sh

This will run the program once and store all intermediate images from the algorithm.
