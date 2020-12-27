# FRF_AMfe_gmsh
Frequency response from mesh generated from Gmsh using the AMfe python library

# Dependencies

   - `numpy`, `scipy`, `pandas`, `matplotlib`
   - `AMfe`, `Gmsh`

# Installation

`AMfe`

A separate environment is recommended in anaconda for your amfe installation. 

Use

git clone https://github.com/AppliedMechanics/AMfe.git

to get the package. For installation of the package in development mode run

cd AMfe
conda create --name <environment-name-of-choice> python=3.7
conda activate <environment-name-of-choice> 
python conda_setup.py
python setup.py develop [no_fortran]

in the main folder. The conda_setup.py file installs the dependencies via conda. It is recommended to install the dependencies with conda because setup.py can only install them via pip which can lead to an unclean conda environment.

The python setup.py develop command builds the fortran routines and installs the python module in-place, i.e., when you do changes to the source code they will be used the next time the module is loaded.

If you do not want to install the FORTRAN-routines, you can add the flag no_fortran to your installation command:

python setup_develop.py develop no_fortran

If no FORTRAN-compiler is found, the installation will work only with the no_fortran-flag.

`Gmsh`

For gmsh, run:

pip install gmsh

For furthur details refer https://pypi.org/project/gmsh/.

# References

A delayed frequency preconditioner approach for speeding-up frequency response computatio of structural components.
Eccomas Proceedia ID: 9005 / Conference Proceeding ID: 19155 / DOI: 10.47964/1120.9005.19155
Authors: Guilherme Jenovencio, Arul Sivasankar, Zeeshan Saeed, Daniel Rixen 
