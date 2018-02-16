Time-scale phase-weighted stack
===============================

[![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.1154588.svg)](https://doi.org/10.5281/zenodo.1154588)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

Software to compute the time-scale phase-weighted stack (ts-PWS) using frames
of continuous wavelets, including the two-stage stack and the unbiased
phase coherence strategies. Ventosa et al. (GJI, 2017)

Main features
-------------
 * Implement the ts-PWS using the Morlet wavelet or the Mexican hat wavelets  
   implemented using frames of wavelets for the time-frequency expansion.
 * Implement the two-stage stacking method and the unbiased phase coherence.
 * All parameters are automatically set when w0, Q or number of cycles is  
   used the parameters describing the frame are changed accordingly.

Compilation
-----------
To compile execute "make" in the src directory. Use "make clean" to remove 
any previouly compiled code.

 * The Seismic Analysis Code (SAC) is used to read and write sac files.
 * The SACDIR enviorment variable should provide the path to directory where
   sac is installed. For example, in bash this can be defined as:  
   export SACDIR=/opt/sac  
   Change "/opt/sac" to your current sac directory if necessary.
 * When available, OpenMP is used to speed up computations.

Warming up
----------
 1. Read ./examples/example.sh
 2. Execute it, e.g., bash example.sh
 3. Do ts-pws for the parameters usage.
 
Paper to be cited
-----------------
Ventosa S., Schimmel M. & E. Stutzmann, 2017. Extracting surface waves,
hum and normal modes: Time-scale phase-weighted stack and beyond,
Geophysical Journal International, 211, 30-44, doi: 10.1093/gji/ggx284

Origin of PWS
-------------
Schimmel M. & J. Gallart, 2007. Frequency-dependent phase coherence for 
noise suppression in seismic array data, Journal of Geophysical Research, 
112, B04303, doi: 10.1029/2006JB004680

Schimmel M. & H. Paulssen, 1997. Noise reduction and detection of weak, 
coherent signals through phase weighted stacks, Geophysical Journal 
International, 130, 497-505, doi: 10.1111/j.1365-246X.1997.tb05664.x

2018/02/15 Sergi Ventosa Rahuet (sergiventosa(at)hotmail.com)
