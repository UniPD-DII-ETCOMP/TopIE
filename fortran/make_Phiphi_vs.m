clear 
close all
clc

mex -O -largeArrayDims -output funPphiphi3_vs_for P_funLphiphi3_vs_mex.f90 funPphiphi3_vs.f90 
