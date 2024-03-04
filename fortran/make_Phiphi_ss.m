clear 
close all
clc

mex -O -largeArrayDims -output funPphiphi3_ss_for funPphiphi3_ss_mex.f90 funPphiphi3_ss.f90 
