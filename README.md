# AMRTransport
Adaptive Mesh Refinement (AMR) for Single Tracer Transport Component
The project includes flux-form semi-Lagrangian scheme for passive tracer transport using adaptive mesh refinement.  
This repository includes codes both for idealized tests and ploting and a patch for the ECHAM6-HAMMOZ.  
The current version supports the use in Linux environment with bash shells.

# Prerequisite for idealized test
* Fortran: gfortran version >= 7.0
* HDF5
* SZIP
* NetCDF4
* NetCDF4-Fortran
* Python package: 
  * Pyngl
  * Nio
  * Numpy
  * Matplotlib.pyplot

# Running
  * To generate idealized test results and plots, type following commands in the terminal:  
    `git clone https://github.com/yumengch/AMRTransport.git`  
    `cd AMRTransport`  
    Open the file `run.sh` and edit the LibPath and LibPathReverse to specify the path to the library  
    e.g.:  
    Substitute `(Plase type the path to the include directory of your netcf-fortran. e.g.: \/usr\/include)`
    into
    `\/usr\/include`   
    `bash run.sh`
  * For the realistic test (to be done):
    1. Obtain the ECHAM6-HAMMOZ
    2. Patch `src/realistic/src.patch` to the src directory in the ECHAM6-HAMMOZ
    3. Follow instructions to recompile the ECHAM6-HAMMOZ
    4. Setup-files for experiments can be asked upon request
    5. corresponding plotting routines are in the `src/realistic/plot`


#### If you have any question regarding the project, E-mail: yumeng.chen@uni-hamburg.de