#!/bin/bash

cd src/idealized/
echo $(pwd)
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
