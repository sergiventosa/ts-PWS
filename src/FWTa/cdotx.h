#ifndef CDOTX_H
#define CDOTX_H

#include <complex.h>

void cdotx_dd (double *y, double *x0, unsigned int N, double *w0, unsigned int L, int c, int inc);
void cdotx_dc (double complex *y, double *x0, unsigned int N, double complex *w0, unsigned int L, int c, int inc);
void cdotx_cd (double complex *y, double complex *x0, unsigned int N, double *w0, unsigned int L, int c, int inc);
void cdotx_cc (double complex *y, double complex *x0, unsigned int N, double complex *w0, unsigned int L, int c, int inc);
void re_cdotx_cc (double *y, double complex *x0, unsigned int N, double complex *w0, unsigned int L, int c, int inc);

void cdotx_upsampling_dd (double *y, unsigned int N, double *x0, unsigned int SX, double *w0, unsigned int L, int c, int D);
void cdotx_upsampling_dc (double complex *y, unsigned int N, double *x0, unsigned int SX, double complex *w0, unsigned int L, int c, int D);
void cdotx_upsampling_cd (double complex *y, unsigned int N, double complex *x0, unsigned int SX, double *w0, unsigned int L, int c, int D);
void cdotx_upsampling_cc (double complex *y, unsigned int N, double complex *x0, unsigned int SX, double complex *w0, unsigned int L, int c, int D);
void re_cdotx_upsampling_cc (double *y, const unsigned int N, double complex *x0, unsigned int SX, const double complex *w, const unsigned int L, const int c, const int D);

#endif
