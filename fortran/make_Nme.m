clear 
close all
clc
mex -O -largeArrayDims -output funNme_for funNme_mex.f90 funNme.f90 
