#ifndef MATLABIOWT_H
#define MATLABIOWT_H

#ifdef MATLAB
#include <complex.h>
#include "mex.h"
#include "wavelet_v7.h"

mxArray *export_RealWaveletDec(double **d, unsigned int *N, unsigned int S);
mxArray *export_ComplexWaveletDec(double complex **d, unsigned int *N, unsigned int S);
mxArray *export_RWT(t_RWTvar *x, int M, int N);
mxArray *export_CWT(t_CWTvar *x, int M, int N);
mxArray *export_WaveletFamily(t_WaveletFamily *pWF);
mxArray *export_WaveletFamilyArray(t_WaveletFamily **pWF, int M, int N);

#endif

#endif
