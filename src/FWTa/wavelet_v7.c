/*
 * Minimal wavelet library using frames. Based on the wavelet library developed in my PhD. for the time-scale slowness filters.
 * - Main decompostion and reconstruction functions.
 * 2011/01 Sergi Ventosa
 * 2016/01 Sergi Ventosa: migrated to c99 double complex types.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "wavelet_v7.h"
#include "cdotx.h"
#include "myallocs.h"
#include "prnmsg.h"

/********************************************************************************/
/* Continuous Wavelet Transform for real and complex signals and their inverses */
/* including optional down/upsamplings.                                         */
/********************************************************************************/
int real_1D_wavelet_dec (t_RWTvar *y, double *x, unsigned int N, t_WaveletFamily *pWF) {
	unsigned int s, L;
	int er = 0;

	if (x == NULL) {
		prerror("real_1D_wavelet_dec: data null pointer.");
		return -1;
	}
	y->S = pWF->Ns;
	#pragma omp parallel for default(shared) private(L) schedule(dynamic,1)
	for (s=0; s<pWF->Ns; s++) {
		L = pWF->Ls[s];
		if (N < L) {
			prerror("real_1D_wavelet_dec: scale = %d, Ls = %d > N = %d", s, pWF->Ls[s], N);
			er = -3; 
			continue;
		}
		cdotx_dd(y->d[s], x, N, pWF->wframe.wr[s], L, pWF->center[s], pWF->Down_smp[s]);
		y->N[s] = (N + pWF->Down_smp[s] - 1)/pWF->Down_smp[s];
	}
	return er;
}

int complex_1D_wavelet_dec (t_CWTvar *y, double *x, unsigned int N, t_WaveletFamily *pWF) {
	unsigned int s, L;
	int er = 0;

	if (x == NULL) {
		prerror("complex_1D_wavelet_dec: data null pointer.");
		return -1;
	}
	y->S = pWF->Ns;
	// #pragma omp parallel for default(shared) private(L) schedule(dynamic,1)
	for (s=0; s<pWF->Ns; s++) {
		L = pWF->Ls[s];
		if (N < L) {
			prerror("complex_1D_wavelet_dec: scale = %d, Ls = %d > N = %d", s, pWF->Ls[s], N);
			er = -3; 
			continue;
		}
		cdotx_dc(y->d[s], x, N, pWF->wframe.wc[s], L, pWF->center[s], pWF->Down_smp[s]);
		y->N[s] = (N + pWF->Down_smp[s] - 1)/pWF->Down_smp[s];
	}
	return er;
}

int real_1D_wavelet_rec (double *xrec, t_RWTvar *y, unsigned int N, t_WaveletFamily *pWF) {
	double *bf0, *bf, *p, aux;
	unsigned int u, s;

	if (y == NULL) {
		prerror("real_1D_wavelet_rec: RWTvar null pointer.");
		return -1;
	}
	p = xrec;
	for (u=0; u<N; u++) *p++ = 0;

	if (NULL == (bf0 = (double *)mymalloc(N*sizeof(double)))) {
		prerror("real_1D_wavelet_rec: run out of memory.");
		return -2;
	}
	for (s=0; s<pWF->Ns; s++) {
		bf = bf0;
		p = xrec;
		if (pWF->Down_smp[s] > 1)
			cdotx_upsampling_dd(bf, N, y->d[s], y->N[s], pWF->wdualframe.wr[s], pWF->Lds[s], pWF->center_df[s], pWF->Down_smp[s]);
		else 
			cdotx_dd(bf, y->d[s], N, pWF->wdualframe.wr[s], pWF->Lds[s], pWF->center_df[s], 1);
		aux = log(pWF->a0) / (pWF->Cpsi * pWF->V * pWF->scale[s]);
		for (u=0; u<N; u++) p[u] += aux * bf[u];
	}
	myfree(bf0);
	return 0;
}

int complex_1D_wavelet_rec (double complex *xrec, t_CWTvar *y, unsigned int N, t_WaveletFamily *pWF) {
	double complex *p, *bf0, *bf;
	double aux;
	unsigned int u, s;

	if (y == NULL) {
		prerror("complex_1D_wavelet_rec: CWTvar null pointer.");
		return -1;
	}
	for (u=0; u<N; u++) xrec[u] = 0;
	
	if (NULL == (bf0 = (double complex *)mymalloc(N*sizeof(double complex)) )) {
		prerror("complex_1D_wavelet_rec: run out of memory.");
		return -2;
	}
	for (s=0; s<pWF->Ns; s++) {
		bf = bf0;
		p = xrec;
		if (pWF->Down_smp[s] > 1)
			cdotx_upsampling_cc(bf, N, y->d[s], y->N[s], pWF->wdualframe.wc[s], pWF->Lds[s], pWF->center_df[s], pWF->Down_smp[s]);
		else
			cdotx_cc(bf, y->d[s], N, pWF->wdualframe.wc[s], pWF->Lds[s], pWF->center_df[s], 1);
		aux = log(pWF->a0) / (2 * pWF->Cpsi * pWF->V * pWF->scale[s]);
		for (u=0; u<N; u++)	p[u] += aux * bf[u];
	}
	myfree(bf0);
	return 0;
}

int Re_complex_1D_wavelet_rec (double *xrec, t_CWTvar *y, unsigned int N, t_WaveletFamily *pWF) {
	double *bf0, *bf, *p, aux;
	unsigned int u, s;

	if (y == NULL) {
		prerror("re_complex_1D_wavelet_rec: CWTvar null pointer.");
		return -1;
	}
	for (u=0; u<N; u++) xrec[u] = 0;

	if (NULL == (bf0 = (double *)mymalloc(N*sizeof(double)) )) {
		prerror("Re_complex_1D_wavelet_rec: run out of memory.");
		return -2;
	}
	for (s=0; s<pWF->Ns; s++) {
		bf = bf0;
		p = xrec;
		if (pWF->Down_smp[s] > 1)
			re_cdotx_upsampling_cc(bf, N, y->d[s], y->N[s], pWF->wdualframe.wc[s], pWF->Lds[s], pWF->center_df[s], pWF->Down_smp[s]);
		else
			re_cdotx_cc(bf, y->d[s], N, pWF->wdualframe.wc[s], pWF->Lds[s], pWF->center_df[s], 1);
		aux = log(pWF->a0) / (2 * pWF->Cpsi * pWF->V * pWF->scale[s]);
		for (u=0; u<N; u++) p[u] += aux * bf[u];
	}
	myfree(bf0);
	return 0;
}

int Nreal_1D_wavelet_dec(t_RWTvar *y, double *x, unsigned int N, unsigned int M, t_WaveletFamily *pWF) {
	unsigned int i;
	int er=0;

	if (x == NULL) {
		prerror("Nreal_1D_wavelet_dec: data null pointer.");
		return -1;
	}
	for (i=0; i < M; i++) {
		er = real_1D_wavelet_dec (y++, x, N, pWF);
        if (er) {
			prerror("Nreal_1D_wavelet_dec: step %d of %d", i, M);
			break;
		}
		x += N;
	}
    return er;
}

int Nreal_1D_wavelet_rec(double *xrec, t_RWTvar *y, unsigned int N, unsigned int M, t_WaveletFamily *pWF) {
	unsigned int i;
	int er=0;

	if (y == NULL) {
		prerror("Ncomplex_1D_wavelet_rec: RWTvar null pointer.");
		return -1;
	}
	for (i=0; i < M; i++) {
		er = real_1D_wavelet_rec (xrec, y++, N, pWF);
		if (er) {
			prerror("Ncomplex_1D_wavelet_rec: step %d of %d", i, M);
			break;
		}
		xrec += N;
	}
    return er;
}

int Ncomplex_1D_wavelet_dec(t_CWTvar *y, double *x, unsigned int N, unsigned int M, t_WaveletFamily *pWF) {
	unsigned int i;
	int er=0;

	if (x == NULL) {
		prerror("Ncomplex_1D_wavelet_dec: data null pointer.");
		return -1;
	}
	for (i=0; i < M; i++) {
		er = complex_1D_wavelet_dec (y++, x, N, pWF);
        if (er) {
			prerror("Ncomplex_1D_wavelet_dec: step %d of %d", i, M);
			break;
		}
		x += N;
	}
    return er;
}

int Ncomplex_1D_wavelet_rec(double complex *xrec, t_CWTvar *y, unsigned int N, unsigned int M, t_WaveletFamily *pWF) {
	unsigned int i;
	int er=0;

	if (y == NULL) {
		prerror("Ncomplex_1D_wavelet_rec: CWTvar null pointer.");
		return -1;
	}
	for (i=0; i < M; i++) {
		er = complex_1D_wavelet_rec (xrec, y++, N, pWF);
		if (er) {
			prerror("Ncomplex_1D_wavelet_rec: step %d of %d", i, M);
			break;
		}
		xrec += N;
	}
    return er;
}

int NRe_complex_1D_wavelet_rec(double *xrec, t_CWTvar *y, unsigned int N, unsigned int M, t_WaveletFamily *pWF) {
	unsigned int i;
	int er=0;

	if (y == NULL) {
		prerror("Nre_complex_1D_wavelet_rec: CWTvar null pointer.");
		return -1;
	}
	for (i=0; i < M; i++) {
		er = Re_complex_1D_wavelet_rec (xrec, y++, N, pWF);
		if (er) {
			prerror("Nre_complex_1D_wavelet_rec: step %d of %d", i, M);
			break;
		}
		xrec += N;
	}
    return er;
}

