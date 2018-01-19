Time-scale phase-weighted stack
===============================

[![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.1154588.svg)](https://doi.org/10.5281/zenodo.1154588)

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
 
Papers to be cited
------------------
Ventosa, S., Schimmel, M., & Stutzmann, E., 2017. Extracting surface waves,
hum and normal modes: Time-scale phase-weighted stack and beyond,
Geophysical Journal International, 211, 30-44, doi: 10.1093/gji/ggx284

2017/09/20 Sergi Ventosa Rahuet (sergiventosa(at)hotmail.com)
