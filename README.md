<p align="center">
	<img src="image.png" width="486.4">
</p>

# TopIE 

This directory contains the TopIE code for integral-equation-based topology optimization

-------------------------------------------------------------------

# Description
 
TopIE.m is the main file you must run to start the code. 

All user-settable quantities, e.g. current density, frequency, TOBS settings are contained in the block identified by the 
BEGIN USER SETTINGS / END USER SETTINGS comments.

Available test cases
--------------------
Some simple test cases are contained in separate directories under "test_cases". 
Set the "test_case_dir" variables in "TopIE.m"  to the appropriate directories.

User-defined test cases
-----------------------
Follow the instuctions given in "README.txt" inside the "test_cases" directory.

Results visualization
--------------------
The plot of the power losses is given at the end of the simulation.

Credits
--------------------
If you use TopIE, please consider citing:

 [1] [F. Lucchini et al., "TopIE: An Integral Equation Tool for Topology Optimization in Electromagnetics," in IEEE Transactions on Antennas and Propagation, vol. 72, no. 1, pp. 1-8, Jan. 2024, doi: 10.1109/TAP.2023.3321143](https://ieeexplore.ieee.org/document/10273800)
 
and TopIE itself

 [2] F. Lucchini, R. Torchio, "TopIE toolbox", https://github.com/UniPD-DII-ETCOMP/TopIE

Contacts & Authors
-----------------------
For any questions or need help to solve your problem, please contact us

Francesco Lucchini (francesco.lucchini@unipd.it)

Riccardo Torchio (riccardo.torchio@unipd.it)
