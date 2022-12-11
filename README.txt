Guidelines how to use the code: 

REQUIREMENTS:
 - Cmake version no less than 3.7.2
 - A suitable C++ compiler
 - The linear algebra library Eigen
 - The Boost libraries
 
 
Further the code uses the open source C++-package on cubic spline interpolation included in the class spline.h which is from:
https://kluge.in-chemnitz.de/opensource/spline/ 


After compilation using cmake and make, a target called "s2d" is produced 
which can be executed by: ode 

By this the 1d-simulation is started for the surrogate model for a cylinder domain O_h

The main code of the simulation is provided in the executable 
file "ode.cpp".

In "ode.cpp" the relevant parameters are described 
in lines 34-41, where e.g. the different parameters for the elastic moduli
and the boundary conditions are set

The default configuration of the code is appropriate to
the 1d-simulation provided in the thesis Section 15.
