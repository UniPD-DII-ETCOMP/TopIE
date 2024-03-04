clear 
close all
clc
disp('Do it out from dropbox')

mex -O -largeArrayDims -output gcd_mexed_HEXA commonstuff_hexa.f90 gcd_hexa_mex.f90 gcd_hexa.f90 meshtopoaux.f90

