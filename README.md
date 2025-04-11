FORTRAN code for inverting co-seismic surface deformation data (GPS, InSAR, etc.) for fault slip distribution.

Highlights:

(1) incorporating with layered crust structure

(2) arbitrarily curved fault geometry

(3) a-priori constraint on the variation range of rake angle

(4) optional smoothing constraint applied to slip or stress-drop

(5) fast optimization algorithm based on the steepest descent method 

For Windows user, the executable file is provided under the folder WindowsEXE. Linux user may compile the source codes with "gfortran" via a single command like, e.g.,
~>cd .../SDM2011/SourceCode
~>gfortran -o sdm2011 *.f -O3
to get the excutable code sdm2011.
