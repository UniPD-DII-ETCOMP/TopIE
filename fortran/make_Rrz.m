clear 
close all
clc
mex -O -largeArrayDims -output funRrz_for funRrz_mex.f90 funRrz.f90 
