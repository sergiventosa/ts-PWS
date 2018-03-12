/*****************************************************************************/
/* Program to perform time-scale phase-weiged stack (ts-pws) using frames    */
/* of continuous wavelets.                                                   */
/*                                                                           */
/* Author: Sergi Ventosa Rahuet (sergiventosa@hotmail.com)                   */
/*                                                                           */
/* The actual library computing ts_pws & friends.                            */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "wavelet_v7.h"
#include "myallocs.h"
#include "prnmsg.h"
#include "ts_pws1e_lib.h"

double similarity (double *x1, float *x2, unsigned int N);
double misfit (double *x1, float *x2, unsigned int N);
int tspws_stacks_float (t_CWTvar *pWTst, t_CWTvar *pWTps, float *sigall, int max, 
	unsigned int nsamp, unsigned int mtr, t_WaveletFamily *pWF, t_CWTvar *pWT);
int tspws_stacks_float_1step (t_CWTvar *pWTst, t_CWTvar *pWTps, float *sigall, int max, 
	unsigned int nsamp, unsigned int itr, t_WaveletFamily *pWF, t_CWTvar *pWT);
void partial_linear_stacks (double *dsigN, float *sigall, unsigned int max, unsigned int nsamp, 
	unsigned int Kmax, unsigned int mtr);
void tspws_stacks_double (t_CWTvar *pWTst, t_CWTvar *pWTps, double *dsigall, int max, 
	unsigned int nsamp, unsigned int mtr, t_WaveletFamily *pWF, t_CWTvar *pWT);
void tspws_biased (t_CWTvar *pWTout, t_CWTvar *pWTst, t_CWTvar *pWTps, unsigned int nsamp, 
	unsigned int Kmax, unsigned int mtr, double wu);
void tspws_unbiased (t_CWTvar *pWTout, t_CWTvar *pWTst, t_CWTvar *pWTps, unsigned int msamp, 
	unsigned int Kmax, unsigned int mtr);

int tspws_main(t_tsPWS *tspws, t_tsPWS_out *out, t_data *in) {
	float *sigall, *sig, *sig2;
	float fa1, beg, dt;
	double rmean, da1, *dsig;
	unsigned int n, nsmp, nsamp, itr, mtr;
	int nerr=0, max;
	t_hdr *hdr;
	t_WaveletFamily *pWF=NULL;
	
	if (tspws == NULL || out == NULL || in == NULL ) { printf("tspws_main: NULL input\n"); return -1; }

	/********************/
	/* Initializations. */
	/********************/
	sigall = in->sigall;
	hdr    = &in->hdr;
	max    = hdr->max;
	mtr    = (tspws->Nmax) ? tspws->Nmax : hdr->mtr;
	beg    = hdr->beg;
	dt     = hdr->dt;
	nsamp  = max;
	
	/* Folding */
	if (tspws->fold) { /* Sequences are not halved to avoid the jump between 0 and -beg in the WT. */
		if (2*beg + (max-1)*dt > 0.5*dt) {
			printf("Warning: Folding ignored. B = %f, E = %f, nsamp = %u\n", beg, beg + (max-1)*dt, nsamp);
			tspws->fold = 0;
		} else {
			nsmp = nsamp/2;  /* Interger (floor) division. */
			for (itr=0; itr<mtr; itr++) {
				sig = sigall + itr*max;
				for (n=0; n<nsmp; n++) {
					fa1  = sig[n];  
					fa1 += sig[max-1-n];
					fa1 *= 0.5; 
					sig[max-1-n] = fa1;
				}
			} 
		}
	}
	
	/* Wavelet initialization */
	if (tspws->fmin != 0 && tspws->fmin < 1/(dt*nsamp)) {
		printf("Warning: fmin is too low. Replaced by the default value.\n");
		tspws->fmin = 0;
	}
	switch (tspws->w0set) { /* Valid for the Morlet wavelet, ignored when using mexican hat. */
		case 1: tspws->w0 = 2*sqrt(log(2))  * tspws->Q;     break;
		case 2: tspws->w0 = PI/sqrt(log(2)) * tspws->cycle; break;
	}
	if (tspws->type == -1 || tspws->type == -2) {
		da1 = tspws->w0 / (PI*sqrt(2/log(2)));
		if (!tspws->lVfix)  tspws->V  = (unsigned int)ceil( 4. * da1 );
		if (!tspws->lb0fix) tspws->b0 = (unsigned int)pow( 2, round(log2(da1)) );
		if (!tspws->ls0fix) tspws->s0 = 2.;
	}
	if (tspws->type == -3) {
		tspws->w0 = sqrt(2);
		if (!tspws->lVfix)  tspws->V  = 2;
		if (!tspws->lb0fix) tspws->b0 = 0.5;
		if (!tspws->ls0fix) tspws->s0 = 1.;
	}
	if (tspws->fmin) {
		da1 = tspws->w0/(2*PI * dt*tspws->fmin);  /* Fmin sets the largest scale.       */
		if (tspws->J) {                           /* J is defined:                      */
			da1 /= pow(2, tspws->J - 1/(double)tspws->V); /*   da1 is the now the lowest scale. */
			while (da1 < tspws->s0 * 0.9) {             /*   Check that it's not too low.     */
				da1 *= 2;
				tspws->J--;
			}
			tspws->s0 = da1;
		} else tspws->J = (unsigned int)floor(log2(da1/tspws->s0) + 1/(double)tspws->V);
	} else if (!tspws->J) {
		da1 = nsamp*tspws->w0 / (2*PI * 4.*tspws->s0);
		tspws->J = (unsigned int)floor(log2(da1) + 1/(double)tspws->V);
	}
	
	if (tspws->verbose) {
		printf("Sequence length = %d, dt = %f\n", nsamp, dt);
		printf("  PWS power = %f\n", tspws->wu);
		printf("  Two-stage ");
		if(tspws->Kmax) printf("on, %d groups\n", tspws->Kmax);
		else printf("off\n");
		printf("  unbiased "); if(tspws->unbiased) printf("on\n"); else printf("off\n");
		printf("  rm ");       if(tspws->lrm)      printf("on\n"); else printf("off\n");
		printf("  fold ");     if(tspws->fold)     printf("on\n"); else printf("off\n");
		printf("Mother wavelet: ");
		switch (tspws->type) {
			case -1: printf("complex Morlet\n");       break;
			case -2: printf("exact complex Morlet\n"); break;
			case -3: printf("complex Mexican hat\n");  break;
		}
		printf("Sampling of the time-frequency domain:\n");
		printf("         V = %d, J = %d, b0 = %f, s0 = %f, w0 = %f\n", tspws->V, tspws->J, tspws->b0, tspws->s0, tspws->w0);
		printf("This is: fmax = %f Hz, fmin = %f Hz\n", tspws->w0/(2*PI*tspws->s0*dt), 
			tspws->w0/(2*PI*dt*tspws->s0 * pow(2, tspws->J - 1/(double)tspws->V)) );
		printf("         Q = %f (equivalently cycles = %f)\n\n", tspws->w0/(2*sqrt(log(2))), tspws->w0*sqrt(log(2))/PI);
	}
	
	pWF = CreateWaveletFamily (tspws->type, tspws->J, tspws->V, nsamp, tspws->s0, tspws->b0, 0, tspws->w0, tspws->uni);
	
	/*******************/
	/* Main Processing */
	/*******************/
	/* remove the mean */
	#pragma omp parallel for default(shared) private(sig, n, fa1, rmean) schedule(static)
	for (itr=0; itr<mtr; itr++) {
		sig = sigall + itr*max;
		if (tspws->lrm) {
			rmean = 0.;
			/* Doubles are needed because nsamp can be large and then high numerical errors can be made. */
			for(n=0; n<nsamp; n++) rmean += (double)sig[n];
			fa1 = (float)(rmean/nsamp);
			for(n=0; n<nsamp; n++) sig[n] -= fa1;
		}
	}
	
	/*****************************/
	/* The linear t-domain stack */
	/*****************************/
	sig2 = out->ls;
	for (itr=0; itr<mtr; itr++) {
		sig = sigall + itr*max;
		for(n=0; n<nsamp; n++) sig2[n] += sig[n];
	}
	fa1 = 1.0/mtr;
	for(n=0; n<nsamp; n++) sig2[n] *= fa1;
	
	/**************************************/
	/* The time-scale phase weighed stack */
	/**************************************/
	if (mtr) {
		t_CWTvar *pWT, *pWTst, *pWTps, *pWTout;
		unsigned int Kmax;
		double *dsigN;
		
		if (NULL == (pWT    = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The trace in time-scale         */
		if (NULL == (pWTst  = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The stacked part (& the result) */
		if (NULL == (pWTps  = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The phase weight part.          */
		if (NULL == (pWTout = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The phase weight part.          */
		if (NULL == (dsig   = (double *)mymalloc(max * sizeof(double)) )) nerr = 4; /* The ts-pws result.              */
	
		/* Linear and phase stacks. */
		Kmax = tspws->Kmax;
		if (!nerr) {
			if (!tspws->Kmax || tspws->Kmax > mtr) {
				Kmax = mtr;
				tspws_stacks_float (pWTst, pWTps, sigall, max, nsamp, mtr, pWF, pWT);
			} else {
				if (NULL == (dsigN = (double *)malloc(Kmax*nsamp * sizeof(double)) )) nerr = 4;
				else {
					/* Kmax partials linear stacks */
					partial_linear_stacks (dsigN, sigall, max, nsamp, Kmax, mtr);
					/* Linear and phase stack of the partial stacks. */
					tspws_stacks_double (pWTst, pWTps, dsigN, max, nsamp, Kmax, pWF, pWT);
					/* Some cleaning */
					myfree(dsigN);
				}
			}
			
			/* Out of the loop part of the phase-weighed stack. */
			/* The actual ts-PWS. */
			if (tspws->wu == 2 && tspws->unbiased) 
				tspws_unbiased (pWTout, pWTst, pWTps, nsamp, Kmax, mtr);
			else tspws_biased (pWTout, pWTst, pWTps, nsamp, Kmax, mtr, tspws->wu);
			/* ICWT */
			Re_complex_1D_wavelet_rec (dsig, pWTout, nsamp, pWF);
			/* Save the results */
			sig2 = out->tsPWS;
			for (n=0; n<nsamp; n++) sig2[n] = (float)dsig[n];
			for (   ; n<max;   n++) sig2[n] = 0.;
			myfree(dsig);
		}
		
		/***************/
		/* Convergence */
		/***************/
		if (tspws->convergence) {
			unsigned int Tr;
			float *pf1;
			
			if (NULL == ( dsig    = (double *)mymalloc(max * sizeof(double)) )) nerr = 4;
			if (NULL == (dsigN    = (double *)mymalloc(Kmax*max * sizeof(double)) )) nerr = 4;
			
			if (!nerr) {
				/* ts-PWS convergence */
				CleanComplexWaveletVar (pWTst);
				CleanComplexWaveletVar (pWTps);
				pf1 = out->tsPWS_steps;
				for (itr=0; itr<mtr; itr++) {
					Tr = itr + 1;
					/* ts-PWS on itr traces, about the same code as above */
			 		if (!tspws->Kmax || tspws->Kmax >= Tr) {
						Kmax = Tr;
						tspws_stacks_float_1step (pWTst, pWTps, sigall, max, nsamp, itr, pWF, pWT);
					} else {
						Kmax = tspws->Kmax;
						partial_linear_stacks (dsigN, sigall, max, nsamp, Kmax, Tr);
						tspws_stacks_double (pWTst, pWTps, dsigN, max, nsamp, Kmax, pWF, pWT);
					}
					
					if (tspws->wu == 2 && tspws->unbiased) 
						tspws_unbiased (pWTout, pWTst, pWTps, nsamp, Kmax, Tr);
					else tspws_biased (pWTout, pWTst, pWTps, nsamp, Kmax, Tr, tspws->wu);
					/* ICWT */
					Re_complex_1D_wavelet_rec (dsig, pWTout, nsamp, pWF);
					/* Convergence curves */
					sig2 = (in->reference) ? in->reference : out->tsPWS;
					out->tsPWS_sim[itr]    = similarity (dsig, sig2, nsamp);
					out->tsPWS_misfit[itr] = misfit (dsig, sig2, nsamp);
					/* Save the stacked traces for each itr. */
					if (out->tsPWS_steps) {
						for (n=0; n<nsamp; n++) pf1[n] = (float)dsig[n];
						for (   ; n<max;   n++) pf1[n] = 0.;
						pf1 += max;
					}
				}
				
				/* Linear-stack convergence */
				pf1 = out->ls_steps;
				for(n=0; n<nsamp; n++) dsig[n] = 0;
				for (itr=0; itr<mtr; itr++) {
					sig = sigall + itr*max;
					for(n=0; n<nsamp; n++) dsig[n] += sig[n];
					
					fa1 = 1.0/(itr+1);
					for(n=0; n<nsamp; n++) dsig[n] *= fa1;
					/* Convergence curves */
					sig2 = (in->reference) ? in->reference : out->ls;
					out->ls_sim[itr]    = similarity (dsig, sig2, nsamp);
					out->ls_misfit[itr] = misfit (dsig, sig2, nsamp);
					/* Save the stacked traces for each itr. */
					if (out->ls_steps) {
						for (n=0; n<nsamp; n++) pf1[n] = (float)dsig[n];
						for (   ; n<max;   n++) pf1[n] = 0.;
						pf1 += max;
					}
					for(n=0; n<nsamp; n++) dsig[n] *= (itr+1);
				}
			}
			
			/* Some cleaning */
			myfree(dsigN);
			myfree(dsig);
		}
		
		/* Some cleaning */
		DestroyComplexWaveletVar (pWTps);
		DestroyComplexWaveletVar (pWTst);
		DestroyComplexWaveletVar (pWTout);
		DestroyComplexWaveletVar (pWT);
	}
	
	DestroyWaveletFamily (pWF);
	
	return nerr;
}
/* Compute similiarity as GNCC. */
double similarity (double *x1, float *x2, unsigned int N) {
	double y, da1;
	unsigned int n;
	
	y = 0;
	for (n=0; n<N; n++)   y += x1[n] * x2[n];
	
	da1 = 0;
	for (n=0; n<N; n++) da1 += x1[n] * x1[n];
	y /= sqrt(da1);
	
	da1 = 0;
	for (n=0; n<N; n++) da1 += x2[n] * x2[n];
	y /= sqrt(da1);
	
	return y;
}

/* Compute misfit as sum (|x1 - x2|^2). */
double misfit (double *x1, float *x2, unsigned int N) {
	double y, da1;
	unsigned int n;
	
	y = 0;
	for (n=0; n<N; n++) {
		da1 = (x1[n] - (double)x2[n]);
		y += da1*da1;
	}
	return y;
}

/* Compute the linear stack (pWTst) and the phase stack (pWTps) in the time-scale domain of  */
/* the mtr sequences of sigall (float precision).                                            */
int tspws_stacks_float (t_CWTvar *pWTst, t_CWTvar *pWTps, float *sigall, int max, 
			unsigned int nsamp, unsigned int mtr, t_WaveletFamily *pWF, t_CWTvar *pWT) {
	double complex ca1;
	double *dsig;
	unsigned int itr, s, n;
	float *sig;
	
	CleanComplexWaveletVar (pWTst);
	CleanComplexWaveletVar (pWTps);
	dsig = (double *)mymalloc(max * sizeof(double));
	if (dsig == NULL) return 4;
	else {
		for (itr=0; itr<mtr; itr++) {
			sig = sigall + itr*max;
			
			/* WT */
			for(n=0; n<nsamp; n++) dsig[n] = (double)sig[n];
			complex_1D_wavelet_dec (pWT, dsig, nsamp, pWF);
			
			#pragma omp parallel for default(shared) private(n, ca1) schedule(dynamic,1)
			for (s=0; s<pWT->S; s++) {
				for (n=0; n<pWT->N[s]; n++) {
					ca1 = pWT->d[s][n];
					pWTst->d[s][n] += ca1;            /* Linear stack. */
					ca1 /= cabs(ca1);
					if (creal(ca1*conj(ca1)) <= 1.001) pWTps->d[s][n] += ca1;   /* Phase stack.  */
				}
			}
		}
		myfree(dsig);
	}
	return 0;
}

/* Add the itr sequence of sigall to the linear stack (pWTst) and to the phase stack (pWTps) */
/* in the time-scale domain.                                                                 */
int tspws_stacks_float_1step (t_CWTvar *pWTst, t_CWTvar *pWTps, float *sigall, int max, 
			unsigned int nsamp, unsigned int itr, t_WaveletFamily *pWF, t_CWTvar *pWT) {
	double complex ca1;
	double *dsig;
	unsigned int s, n;
	float *sig;
	
	dsig = (double *)mymalloc(max * sizeof(double));
	if (dsig == NULL) return 4;
	else {
		sig = sigall + itr*max;
		
		/* WT */
		for(n=0; n<nsamp; n++) dsig[n] = (double)sig[n];
		complex_1D_wavelet_dec (pWT, dsig, nsamp, pWF);
		
		#pragma omp parallel for default(shared) private(n, ca1) schedule(dynamic,1)
		for (s=0; s<pWT->S; s++) {
			for (n=0; n<pWT->N[s]; n++) {
				ca1 = pWT->d[s][n];
				pWTst->d[s][n] += ca1;            /* Linear stack. */
				ca1 /= cabs(ca1);
				if (creal(ca1*conj(ca1)) <= 1.001) pWTps->d[s][n] += ca1;   /* Phase stack.  */
			}
		}
		myfree(dsig);
	}
	return 0;
}

/* Linearly stack mtr sequences in Kmax groups in the sequencial order. */
void partial_linear_stacks (double *dsigN, float *sigall, unsigned int max, unsigned int nsamp, 
			unsigned int Kmax, unsigned int mtr) {
	unsigned int itr, n;
	float *sig;
	double *pd1;
	
	memset(dsigN, 0, Kmax*max*sizeof(double));
	for (itr=0; itr<mtr; itr++) {
		sig = sigall + itr*max;
		/* pd1 = dsigN + (itr % Kmax)*nsamp; */  /* Interlaced version */
		pd1 = dsigN + (unsigned int)floor((double)(itr*Kmax)/(double)mtr) * max;
		for (n=0; n<nsamp; n++) pd1[n] += (double)sig[n];
	}
}

/* Compute the linear stack (pWTst) and the phase stack (pWTps) in the time-scale domain of  */
/* the mtr sequences of dsigall (double precision).                                           */
void tspws_stacks_double (t_CWTvar *pWTst, t_CWTvar *pWTps, double *dsigall, int max, 
			unsigned int nsamp, unsigned int mtr, t_WaveletFamily *pWF, t_CWTvar *pWT) {
	double complex ca1;
	unsigned int itr, s, n;
	
	CleanComplexWaveletVar (pWTst);
	CleanComplexWaveletVar (pWTps);
	for (itr=0; itr<mtr; itr++) {
		/* WT */
		complex_1D_wavelet_dec (pWT, dsigall + itr*max, nsamp, pWF);
		
		#pragma omp parallel for default(shared) private(n, ca1) schedule(dynamic,1)
		for (s=0; s<pWT->S; s++) {
			for (n=0; n<pWT->N[s]; n++) {
				ca1 = pWT->d[s][n];
				pWTst->d[s][n] += ca1;            /* Linear stack. */
				ca1 /= cabs(ca1);
				if (creal(ca1*conj(ca1)) <= 1.001) pWTps->d[s][n] += ca1;   /* Phase stack.  */
			}
		}
	}
}

/* Compute the PWS in the time-scale domain from the previously computed linear and phase stacks.*/
void tspws_biased (t_CWTvar *pWTout, t_CWTvar *pWTst, t_CWTvar *pWTps, unsigned int nsamp, 
		unsigned int Kmax, unsigned int mtr, double wu) {
	double complex ca1;
	double da1, da2;
	unsigned int n, s;
	
	if (wu == 2) {
		da1 = (double)Kmax;
		da1 = 1/(da1*da1*(double)mtr);
		#pragma omp parallel for default(shared) private(n, ca1, da2) schedule(dynamic,1)
		for (s=0; s<pWTst->S; s++) {
			for (n=0; n<pWTst->N[s]; n++) {
				ca1 = pWTps->d[s][n];
				da2 = (creal(ca1) * creal(ca1) + cimag(ca1) * cimag(ca1)) * da1;
				pWTout->d[s][n] = da2 * pWTst->d[s][n];
			}
		}
	} else {
		#pragma omp parallel for default(shared) private(n, ca1, da2) schedule(dynamic,1)
		for (s=0; s<pWTst->S; s++) {
			for (n=0; n<pWTst->N[s]; n++) {
				ca1 = pWTps->d[s][n];
				da1 = cabs(ca1) / Kmax;
				da1 = pow(da1, wu);
				pWTout->d[s][n] = pWTst->d[s][n] * da1 / mtr;
			}
		}
	}
}

/* Compute the unbiased PWS in the time-scale domain from the previously computed linear and  */
/* phase stacks. The wu parameter is ignored and set to the default value of 2 !!!!           */
void tspws_unbiased (t_CWTvar *pWTout, t_CWTvar *pWTst, t_CWTvar *pWTps, unsigned int nsamp, 
			unsigned int Kmax, unsigned int mtr) {
	double complex ca1;
	double da1;
	unsigned int n, s;
	
	#pragma omp parallel for default(shared) private(n, ca1, da1) schedule(dynamic,1)
	for (s=0; s<pWTst->S; s++) {
		for (n=0; n<pWTst->N[s]; n++) {
			ca1 = pWTps->d[s][n] / Kmax;
			da1 = creal(ca1) * creal(ca1) + cimag(ca1) * cimag(ca1);
			da1 = (Kmax*da1 - 1)/(Kmax - 1);
			pWTout->d[s][n] = pWTst->d[s][n] * da1 / mtr;
		}
	}
}
