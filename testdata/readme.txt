# Generation of test data

## Synopsis:

$ python genTestData.py <testConfigFile> <refImpl> <sourcePath>
<refPath>

 <testConfigFile> is a text file containing list of reference image
 files.

 <refImpl> is the directory where the MATLAB-reference implementation
 resides.

 <sourcePath> the path to the directory where the source images are
 located (and possibly, where resized images are stored).

 <refPath> the path to the directory where the generated references
 images will be written to.

Example content of a testConfigFile:
 
 Name-3-1-1-1.tif
 Name_50-3-1-1-1.tif
 Name_50-3-0.5-0.5-0.5.tif

 Each file name must have the following structure:
 
 <name>[_<scaling>]-<M>-<Wc>-<Ws>-<We>.<ext>

The use of getTestData.py is best explained by example: Consider the
second entry in the above list: Name_50-3-1-1-1.tif. When reading this
entry, the script looks for the files 

 <sourcePath>/Name.0.tif
 <sourcePath>/Name.1.tif
 ...
 <sourcePath>/Name.N.tif

If the resized files do not already exist, the script generated the
following files:

 <sourcePath>/Name_50.0.tif
 <sourcePath>/Name_50.1.tif
 ...
 <sourcePath>/Name_50.N.tif

Each of these files is a corresponding resized image (the scaling
factor being 50% in this case).

Next, since M is 3 in this case, the first three of the above files
are used as an input for the exposure fusion algorithm. Hence, the
script executes

 $ octave expFuse.m <sourcePath> <refPath> A_50 3 1 1 1

in the directory <refImpl>.

This generate a file with the file
     
 <refPath>/Name_50-3-1-1-1.tif 

