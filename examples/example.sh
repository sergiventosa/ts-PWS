#!/bin/sh
# Compile
cd ../src
rm -f ../sac2bin ../ts_pws
make
mv sac2bin ts_pws ..
cd ../examples

# Example 1: Stack a set of sac files
#  1) Create a text file (filelist.txt) with the list of sac files to stack.
#  2) Stack them using the default configuration. 
#    "verbose" is used to show the parameters used.

ls -1 ECH.00Z.CAN.00Z_500days/ECH.00Z.CAN.00Z_ccgn_*.sac > filelist.txt
../ts_pws filelist.txt verbose

# Example 2: Stack the traces listed in filelist.txt with PWS power of 2 (wu=2), starting at
# 4 mHz (fmin=0.004) and ending at 3 octaves (J=3) above. wu=2 set the PWS power to 2, rm
# remove the mean and fold adds negative and positive lag times. 
#  * The "osac" parameter modifies the default name of the output files by adding the text
#    given. Thus, the file containing the linear stacked result would be "tl_example2.sac" 
#    and the ts-PWS "ts_pws_example2.sac". 
#  * The linear stacked result corresponds to the black trace plotted in Fig 5c of 
#    Ventosa et al (GJI 2017) and the ts-PWS result to the blue trace.

../ts_pws filelist.txt osac="example2" wu=2 rm fold fmin=0.004 J=3 verbose

# Example 3: Stack the traces contained in a single binary file with the two-stage stack. 
# The result should be close to the red trace plotted in Fig 5c of Ventosa et al (GJI 2017).
#  1) Use sac2bin code to copy all the traces listed in filelist.txt to a single binary file 
#     named ECH.00Z.CAN.00Z.bin.
#  2) Stack the traces contained in ECH.00Z.CAN.00Z.bin (note the bin parameter) with the 
#     Two-stage ts-PWS method (TwoStage=10) correcting for the biased of pws for wu=2 (unbiased).
#  3) The linear stacked result corresponds to the black trace plotted in Fig 5c of 
#     Ventosa et al (GJI 2017) and the two-stage stack result to the red trace.

../sac2bin filelist.txt ECH.00Z.CAN.00Z.bin
../ts_pws ECH.00Z.CAN.00Z.bin osac="twostage" wu=2 rm bin fold TwoStage=10 unbiased fmin=0.004 J=3 verbose
