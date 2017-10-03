/*
 * Minimal wavelet library using frames. Based on the wavelet library developed in my PhD. for the time-scale slowness filters.
 * 2011/01 Sergi Ventosa
 * 2016/01 Sergi Ventosa: migrated to c99 double complex types.
 */
#ifndef WAVELET_H
#define WAVELET_H

#include <complex.h>

#ifdef MATLAB
#include "mex.h"
#include "matrix.h"
#endif

#ifndef PI
#define	PI 3.14159265358979328
#endif

#define NSIGMAS        5
#define MAXCMPLX_WAVE  3
#define MAXREAL_WAVE   3

typedef struct {
	int Family;         /* Wavelet Family Id.            */
	unsigned int J;     /* Number of scales 2^1 - 2^J.   */
	unsigned int V;	    /* Number of voices per octave.  */
	unsigned int MWL;   /* Maximum wavelet length.       */
	double s0;          /* First scale.                  */
	double b0;          /* Minimum translation step.     */
	double op1;         /* Optional parameter nº1.       */
	int convtype;       /* Circular 0 / Mirror simetric. */
	int continuous;	    /* Continuous (1) / diadic (0).  */
} t_WaveletDef;

typedef struct {
	int format;              /* Real 1 / Complex -1.          */
	int type;                /* Wavelet Family Name Id.       */
	int convtype;            /* Circular 0 / Mirror simetric. */
	union {
		double **wr;
		double complex **wc;
	} wframe, wdualframe;
	double *scale;           /* Vector of scales.             */
	unsigned int *Ls;        /* The length of each wavelet.   */
	unsigned int *Lds;       /* The length of each dual wavelet. */
	int *center;             /* The wavelet origin.           */
	int *center_df;          /* The dual wavelet origin.      */
	unsigned int *Down_smp;  /* Downsampling of each wavelet. */
	unsigned int Ns;         /* Length of vector of scales.   */
	unsigned int V;          /* Voices per octave.            */
	double Cpsi;             /* Normalization constant.       */
	double a0;
	double b0;
	double op1;              /* Only when it have sense.      */
} t_WaveletFamily;

typedef struct {
	double **d;
	unsigned int *N;
	unsigned int S;
} t_RWTvar;

typedef struct {
	double complex **d;
	unsigned int *N;
	unsigned int S;
} t_CWTvar;

t_WaveletFamily *CreateWaveletFamily (int type, unsigned int J, unsigned int V, unsigned int N, double s0, double b0, int convtype, double w0, int continuous);
void DestroyWaveletFamily (t_WaveletFamily *pWF);

t_RWTvar *CreateRealWaveletVar (t_WaveletFamily *pWF, unsigned int N);
t_CWTvar *CreateComplexWaveletVar (t_WaveletFamily *pWF, unsigned int N);
t_RWTvar *CreateRealWaveletVarArray (t_WaveletFamily *pWF, unsigned int N, unsigned int sz);
t_CWTvar *CreateComplexWaveletVarArray (t_WaveletFamily *pWF, unsigned int N, unsigned int sz);
void CleanRealWaveletVar (t_RWTvar *p);
void CleanComplexWaveletVar (t_CWTvar *p);
void CleanRealWaveletVarArray (t_RWTvar *p, unsigned int sz);
void CleanComplexWaveletVarArray (t_CWTvar *p, unsigned int sz);
void DestroyRealWaveletVar (t_RWTvar *p);
void DestroyComplexWaveletVar (t_CWTvar *p);
void DestroyRealWaveletVarArray (t_RWTvar *p0, unsigned int sz);
void DestroyComplexWaveletVarArray (t_CWTvar *p0, unsigned int sz);

int real_1D_wavelet_dec (t_RWTvar *y, double *x, unsigned int N, t_WaveletFamily *pWF);
int complex_1D_wavelet_dec (t_CWTvar *y, double *x, unsigned int N, t_WaveletFamily *pWF);
int real_1D_wavelet_rec (double *xrec, t_RWTvar *y, unsigned int N, t_WaveletFamily *pWF);
int complex_1D_wavelet_rec (double complex *xrec, t_CWTvar *y, unsigned int N, t_WaveletFamily *pWF);
int Re_complex_1D_wavelet_rec (double *xrec, t_CWTvar *y, unsigned int N, t_WaveletFamily *pWF);
int complex_1D_wavelet_rec_test (double complex *xrec, t_CWTvar *y, unsigned int N, t_WaveletFamily *pWF);

int Nreal_1D_wavelet_dec(t_RWTvar *y, double *x, unsigned int N, unsigned int Tr, t_WaveletFamily *pWF);
int Nreal_1D_wavelet_rec(double *xrec, t_RWTvar *y, unsigned int N, unsigned int M, t_WaveletFamily *pWF);
int Ncomplex_1D_wavelet_dec(t_CWTvar *y, double *x, unsigned int N, unsigned int Tr, t_WaveletFamily *pWF);
int Ncomplex_1D_wavelet_rec(double complex *xrec, t_CWTvar *y, unsigned int N, unsigned int M, t_WaveletFamily *pWF);
int NRe_complex_1D_wavelet_rec(double *xrec, t_CWTvar *y, unsigned int N, unsigned int M, t_WaveletFamily *pWF);

#endif

