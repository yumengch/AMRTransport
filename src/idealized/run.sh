#!/bin/sh

gfortran -c -O3 class_kind.f90 class_parameter.f90 class_AMRTypes.f90 class_coord.f90 class_SetupWind.f90 class_face.f90 class_column.f90 class_mesh.f90 class_predict.f90 class_setup.f90 class_refine.f90 class_CISL.f90 class_depart.f90 class_FFSL.f90 class_netcdf.f90 main.f90 -IPathToNetCDFFINC
gfortran -o FFSL_AMR class_kind.o class_parameter.o class_AMRTypes.o class_coord.o class_SetupWind.o class_face.o class_column.o class_mesh.o class_predict.o class_setup.o class_refine.o class_CISL.o \
        class_depart.o class_FFSL.o class_netcdf.o main.o -LPathToNetCDFFLIB -lnetcdff -LPathToSZIPLIB -LPathToHDF5LIB -LPathToNetCDFLIB -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz -lm -ldl

rm *.o
rm *.mod

