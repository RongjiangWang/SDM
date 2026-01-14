FORTRAN code for inverting co-seismic surface deformation data (GPS, InSAR, etc.) for fault slip distribution.

Highlights:

(1) incorporating with Green's functions based on layered crust structure or user's own Green's functions based on 3D crust structure

(2) arbitrarily curved fault geometry

(3) a-priori constraint on the variation range of rake angle

(4) estimate user-defined unknown parameters (offsets in data, systematic errors in data caused by non-seismic sources, e.g., block motion/rotation, spatial linear trends, etc.) simultaneously with the slip inversion

(4) optional smoothing constraint applied to slip or stress-drop

(5) iterative inversion algorithm based on the Steepest Descent Method using the optimized step-size (relaxation) technique provided by Yong Zhang at Peking University (Zhang, person. comm., 2025)

(6) improvement of slip resolution using Zhang's approach (Zhang, 2025)

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../SDM2025/SourceCode

~>gfortran -o sdm2025 *.f -O3

to get the excutable code sdm2025.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

References

Wang, L., R. Wang, F. Roth, B. Enescu, S. Hainzl and S. Ergintav (2009). Afterslip and viscoelastic relaxation following the 1999 M 7.4 Izmit earthquake from GPS measurements, Geophysical Journal International, 178(3), 1220-1237.

Wang, R., B. Schurr, C. Milkereit, Zh. Shao and M. Jin (2011). An improved automatic scheme for empirical baseline correction of digital strong-motion records, Bulletin of the Seismological Society of America, 101(5), 2029–2044.

Wang, R., S. Parolai, M. Ge, M. Jin, T.R. Walter and J. Zschau (2012). The 2011 Mw 9.0 Tohoku Earthquake: Comparison of GPS and Strong-Motion Data. Bulletin of the Seismological Society of America, doi: 10.1785/0120110264. 

Zhang, Y. (2025). A Simple Method for Improving the Resolution of Geodetic Slip Inversion. Geophysical Journal International 241(3), 1781–1790.

------------------------------
Last update: Jan 14, 2026, by Rongjiang Wang
