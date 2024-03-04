clear 
close all
clc
mex -O -largeArrayDims -output funPphiphi3_vv_for P_funLphiphi3_mex.f90 funPphiphi3.f90 
