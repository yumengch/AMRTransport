#!/bin/sh

gfortran -c -O3 class_kind.f90 class_parameter.f90 class_AMRTypes.f90 class_coord.f90 class_SetupWind.f90 class_face.f90 class_column.f90 class_mesh.f90 class_predict.f90 class_setup.f90 class_refine.f90 class_CISL.f90 class_depart.f90 class_FFSL.f90 class_netcdf.f90 main.f90 -I/sw/jessie-x64/netcdf_fortran-4.4.2-static-gcc81/include
gfortran -o FFSL_AMR class_kind.o class_parameter.o class_AMRTypes.o class_coord.o class_SetupWind.o class_face.o class_column.o class_mesh.o class_predict.o class_setup.o class_refine.o class_CISL.o \
        class_depart.o class_FFSL.o class_netcdf.o main.o -L/sw/jessie-x64/netcdf_fortran-4.4.2-static-gcc81/lib -lnetcdff -L/sw/jessie-x64/szip-2.1-static-gccsys/lib -L/sw/jessie-x64/hdf5-1.8.16-static-gccsys/lib -L/sw/jessie-x64/netcdf-4.3.3.1-static-gccsys/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz -lm -ldl

rm *.o
rm *.mod

