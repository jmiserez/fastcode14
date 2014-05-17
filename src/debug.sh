#!/bin/sh
valgrind --track-origin=yes ./bin/driver_01ref_matlab --log file.log --val val.tif 752:50:752 500:50:500 1.0 1.0 1.0 ../testdata/srcImages/A.0.tif ../testdata/srcImages/A.1.tif
