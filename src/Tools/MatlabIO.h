#ifndef MATLABIO_H
#define MATLABIO_H

#ifdef MATLAB
#include <complex.h>
#include "mex.h"

void mexerror(const char *str1, const char *str2);
double complex *DA2CA(const mxArray *mx);
void DA2CA_nomem(double complex *pc, const mxArray *mx);
int CA2DA(mxArray *mx, double complex *xc);
void check_double_scalar(const mxArray *p, char *msg);

#endif

#endif

