Time-scale phase-weighted stack
===============================

[![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.1154588.svg)](https://doi.org/10.5281/zenodo.1154588)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

Software to compute the time-scale phase-weighted stack (ts-PWS), including the 
two-stage stack and the unbiased phase coherence strategies ([Ventosa et al., GJI 2017](https://doi.org/10.1093/gji/ggx284)).

The software packages of [fast phase cross-correlation](https://github.com/sergiventosa/FastPCC) 
and ts-PWS stacking are basic building blocks in the design of efficient 
signal extraction methods from interstation correlations.

Main features
-------------
 * Fast: The time-frequency expansion used in the ts-PWS is implemented using frames
   of wavelets (Morlet and Mexican hat wavelets are available). The most demanding 
   operation is conventionally reading the data.
 * Better filtering: The two-stage stacking and the unbiased phase coherence help
   to reduce signal attenuation and increase noise attenation.
 * Simple setup: All parameters are automatically set, e.g., when w0, Q or number of 
   cycles is changed the parameters describing the frame are adapted accordingly.
 * Resampling techniques: The bootstrap and jackknife techniques can be applied to 
   build a pull of stacked sequences for a subsequent error analysis.

Compilation
-----------
To compile execute "make" in the src directory. Use "make clean" to remove 
any previouly compiled code.

 * The Seismic Analysis Code (SAC) is used to read and write sac files.
 * The SACHOME enviorment variable should provide the path to directory where sac is
   installed. For example, in bash this can be defined as:  
   export SACHOME=/opt/sac  
   Change "/opt/sac" to your current sac directory if necessary.
 * OpenMP is used to speed up computations. When OpenMP is not available, use 
   make -f makefile_NoOpenMP".

Warming up
----------
 1. Read ./examples/example.sh
 2. Execute it, e.g., bash example.sh
 3. Do ts-pws for the parameters usage.
 
Paper to be cited
-----------------
Ventosa S., Schimmel M. & E. Stutzmann, 2017. Extracting surface waves, hum and 
normal modes: Time-scale phase-weighted stack and beyond, Geophysical Journal 
International, 211(1), 30-44, doi: [10.1093/gji/ggx284](https://doi.org/10.1093/gji/ggx284)

Origin of PWS
-------------
Schimmel M. & J. Gallart, 2007. Frequency-dependent phase coherence for noise 
suppression in seismic array data, Journal of Geophysical Research, 112, B04303, 
doi: [10.1029/2006JB004680](https://doi.org/10.1029/2006JB004680)

Schimmel M. & H. Paulssen, 1997. Noise reduction and detection of weak, coherent 
signals through phase weighted stacks, Geophysical Journal International, 130, 
497-505, doi: [10.1111/j.1365-246X.1997.tb05664.x](https://doi.org/10.1111/j.1365-246X.1997.tb05664.x)

2020/10/01 Sergi Ventosa Rahuet (sergiventosa(at)hotmail.com)
