#!/bin/bash
LibPath() {
    sed -i 's/PathToNetCDFFINC/(Plase type the path to the include directory of your netcf-fortran. e.g.: \/usr\/include)/g' $1
    sed -i 's/PathToNetCDFFLIB/(Plase type the path to the lib directory of your netcf-fortran. e.g.: \/usr\/lib)/g' $1
    sed -i 's/PathToSZIPLIB/(Plase type the path to the lib directory of your szip. e.g.: \/usr\/lib)/g' $1
    sed -i 's/PathToHDF5LIB/(Plase type the path to the lib directory of your HDF5. e.g.: \/usr\/lib)/g' $1
    sed -i 's/PathToNetCDFLIB/(Plase type the path to the lib directory of your netcdf. e.g.: \/usr\/lib)/g' $1
}

LibPathReverse() {
    sed -i 's/(Plase type the path to the include directory of your netcf-fortran. e.g.: \/usr\/include)/PathToNetCDFFINC/g' $1
    sed -i 's/(Plase type the path to the lib directory of your netcf-fortran. e.g.: \/usr\/lib)/PathToNetCDFFLIB/g' $1
    sed -i 's/(Plase type the path to the lib directory of your szip. e.g.: \/usr\/lib)/PathToSZIPLIB/g' $1
    sed -i 's/(Plase type the path to the lib directory of your HDF5. e.g.: \/usr\/lib)/PathToHDF5LIB/g' $1
    sed -i 's/(Plase type the path to the lib directory of your netcdf. e.g.: \/usr\/lib)/PathToNetCDFLIB/g' $1
}

LibPath src/idealized/run.sh
LibPath src/idealized/uniform/run.sh
LibPath src/idealized/interp/run.sh
LibPath src/idealized/intermediate/run.sh
python run_script.py
cd plot
python diff.py
python mov_plane.py
python Solid_Eff.py
python test_case_disp.py
python Div_Eff.py
python mov.py
python SolidIllustration.py
python windmov.py
python div_mass.py
python solid_accuracy.py
python SolidTimeSeries.py
LibPathReverse src/idealized/run.sh
LibPathReverse src/idealized/uniform/run.sh
LibPathReverse src/idealized/interp/run.sh
LibPathReverse src/idealized/intermediate/run.sh