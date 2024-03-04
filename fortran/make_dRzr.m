clear 
close all
clc
mex -O -largeArrayDims -output fun_compute_lamdRu_for fun_compute_lamdRu_mex.f90 fun_compute_lamdRu.f90 
