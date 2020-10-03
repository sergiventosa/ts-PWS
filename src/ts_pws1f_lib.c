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
#include "ts_pws1f_lib.h"

double similarity (double *x1, float *x2, unsigned int N);
double misfit (double *x1, float *x2, unsigned int N);
int tspws_stacks_float (t_CWTvar *pWTst, t_CWTvar *pWTps, float *sigall, double *w, int max, 
	unsigned int nsamp, unsigned int mtr, t_WaveletFamily *pWF, t_CWTvar *pWT);
int tspws_subsmpl_float (float **ls_out, float **tsPWS_out, float *sigall, double *w, int max, 
	unsigned int nsamp, unsigned int mtr, unsigned int K, unsigned int M, double wu, int unbiased, 
	t_WaveletFamily *pWF);
int TwoStage_subsmpl_float (float **ls_out, float **tsPWS_out, float *sigall, double *w, int max, 
	unsigned int nsamp, unsigned int mtr, unsigned int Kmax, unsigned int K, unsigned int M, 
	double wu, int unbiased, t_WaveletFamily *pWF);
int tspws_jackknife_float (float **ls_out, float **tsPWS_out, unsigned int *mtr_subsmpl, float *sigall, 
	time_t *time, double *w, int max, unsigned int nsamp, unsigned int mtr, unsigned int K, unsigned int d, 
	unsigned int M, unsigned int C, double wu, int unbiased, t_WaveletFamily *pWF);
int TwoStage_jackknife_float (float **ls_out, float **tsPWS_out, unsigned int *mtr_subsmpl, float *sigall, 
	time_t *time, double *w, int max, unsigned int nsamp, unsigned int mtr, unsigned int Kmax, unsigned int K, 
	unsigned int d, unsigned int M, unsigned int C, double wu, int unbiased, t_WaveletFamily *pWF);
int tspws_stacks_float_1step (t_CWTvar *pWTst, t_CWTvar *pWTps, float *sigall, int max, 
	unsigned int nsamp, unsigned int itr, t_WaveletFamily *pWF, t_CWTvar *pWT);
void partial_linear_stacks (double *dsigN, float *sigall, double *w, unsigned int max, unsigned int nsamp, 
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
	double rmean, da1, *dsig, *w=NULL;
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
					sig[n] = fa1;
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
		printf("         Q = %f (equivalently cycles = %f)\n", tspws->w0/(2*sqrt(log(2))), tspws->w0*sqrt(log(2))/PI);
		if (tspws->subsmpl_N > 0 && tspws->subsmpl_p > 0)
			printf("Subsampling:\n  N = %d, p = %f\n", tspws->subsmpl_N, tspws->subsmpl_p);
		if (tspws->jackknife_n > 0 && tspws->jackknife_d > 0)
			printf("Jackknife Resampling:\n  n = %d, d = %d\n", tspws->jackknife_n, tspws->jackknife_d);
		printf("\n");
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
	/* Done in the time-scale domain for 1) a more fair comparision and 2) error detection.
	sig2 = out->ls;
	if (w == NULL) {
		for (itr=0; itr<mtr; itr++) {
			sig = sigall + itr*max;
			for(n=0; n<nsamp; n++) sig2[n] += sig[n];
		}
		fa1 = 1.0/mtr;
		for(n=0; n<nsamp; n++) sig2[n] *= fa1;
	} else {
		for (itr=0; itr<mtr; itr++) {
			sig = sigall + itr*max;
			fa1 = w[itr]/mtr;
			for(n=0; n<nsamp; n++) sig2[n] += fa1*sig[n];
		}
	}
	*/ 
	/**************************************/
	/* The time-scale phase weighed stack */
	/**************************************/
	if (mtr) {
		t_CWTvar *pWT, *pWTst, *pWTps, *pWTout;
		unsigned int Kmax;
		double *dsigN, *dsig2;
		
		if (NULL == (pWT    = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The trace in time-scale         */
		if (NULL == (pWTst  = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The stacked part (& the result) */
		if (NULL == (pWTps  = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The phase weight part.          */
		if (NULL == (pWTout = CreateComplexWaveletVar (pWF, nsamp) ))     nerr = 4; /* The phase weight part.          */
		if (NULL == (dsig   = (double *)mymalloc(max * sizeof(double)) )) nerr = 4; /* The ts-pws result.              */
		if (NULL == (dsig2  = (double *)mymalloc(max * sizeof(double)) )) nerr = 4; /* The ts-pws result.              */
		
		/* Linear and phase stacks. */
		Kmax = tspws->Kmax;
		if (!nerr) {
			if (!tspws->Kmax || tspws->Kmax > mtr) {
				Kmax = mtr;
				tspws_stacks_float (pWTst, pWTps, sigall, w, max, nsamp, mtr, pWF, pWT);
			} else {
				if (NULL == (dsigN = (double *)malloc(Kmax*nsamp * sizeof(double)) )) nerr = 4;
				else {
					/* Kmax partials linear stacks */
					partial_linear_stacks (dsigN, sigall, w, max, nsamp, Kmax, mtr);
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
			Re_complex_1D_wavelet_rec (dsig2, pWTst, nsamp, pWF);
			/* Save the results */
			sig2 = out->ls;
			for (n=0; n<nsamp; n++) sig2[n] = (float)dsig2[n]/mtr;
			for (   ; n<max;   n++) sig2[n] = 0.;
			myfree(dsig2);
			
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
			
			if (NULL == ( dsig = (double *)mymalloc(max * sizeof(double)) )) nerr = 4;
			if (NULL == (dsigN = (double *)mymalloc(Kmax*max * sizeof(double)) )) nerr = 4;
			
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
						partial_linear_stacks (dsigN, sigall, w, max, nsamp, Kmax, Tr);
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
		
		/***************/
		/* Subsampling */
		/***************/
		if (tspws->subsmpl_N > 0 && tspws->subsmpl_p > 0) {
			unsigned int K = (unsigned int)ceil((double)mtr * tspws->subsmpl_p);
			unsigned int M = tspws->subsmpl_N;
			
			/* Select traces for the submsamples. */
			if (!tspws->Kmax || tspws->Kmax > mtr) /* Single-stage ts-PWS */
				tspws_subsmpl_float (out->ls_subsmpl, out->tsPWS_subsmpl, sigall, w, max, nsamp, mtr, K, M, tspws->wu, tspws->unbiased, pWF);
			else /* Two-stage ts-PWS */
				TwoStage_subsmpl_float (out->ls_subsmpl, out->tsPWS_subsmpl, sigall, w, max, nsamp, mtr, tspws->Kmax, K, M, tspws->wu, tspws->unbiased, pWF);
		}
		
		if (tspws->jackknife_n > 0 && tspws->jackknife_d > 0) {
			unsigned int K = mtr;
			unsigned int d = tspws->jackknife_d;
			unsigned int M = tspws->jackknife_n;
			unsigned int C = out->M;
			
			if (!tspws->Kmax || tspws->Kmax > mtr) /* Single-stage ts-PWS */
				tspws_jackknife_float (out->ls_subsmpl, out->tsPWS_subsmpl, out->mtr_subsmpl, sigall, in->time, w, max, nsamp, mtr, K, d, M, C, tspws->wu, tspws->unbiased, pWF);
			else
				TwoStage_jackknife_float (out->ls_subsmpl, out->tsPWS_subsmpl, out->mtr_subsmpl, sigall, in->time, w, max, nsamp, mtr, tspws->Kmax, K, d, M, C, tspws->wu, tspws->unbiased, pWF);
		}
	}
	
	DestroyWaveletFamily (pWF);
	myfree(w);
	
	return nerr;
}

/* Selects K samples among J */
int SubsamplingPlan (char *sel, unsigned int J, unsigned int K) {
	unsigned int j, k=0;
	
	if (sel == NULL) return 1;
	if (K > J) return 2;
	
	if (2*K < J) { /* More 0s than 1s */
		memset(sel, 0, J*sizeof(char));
		while (k < K) {
			j = rand() % J;
			if (!sel[j]) {
				k++;
				sel[j] = 1;
			}
		}
	} else {      /* More 1s than 0s */
		memset(sel, 1, J*sizeof(char));
		K = J - K;
		while (k < K) {
			j = rand() % J;
			if (sel[j]) {
				k++;
				sel[j] = 0;
			}
		}
	}
	
	return 0;
}

int JackknifePlans (char **selC, time_t *time, unsigned int mtr, unsigned int d, unsigned int n, unsigned int C) {
	unsigned int itr, c, i, j, *setn, **combs, ua1;
	struct tm *tm;
	
	if (selC == NULL || time == NULL) return 1;
	if (time[0] == 0) return -2;
	
	if (NULL == (setn = (unsigned int *)mymalloc(mtr*sizeof(int)) )) return -1;
	if (NULL == (combs = (unsigned int **)mymalloc(C*sizeof(int *)) )) return -1;
	if (NULL == (combs[0] = (unsigned int *)mymalloc(d*C*sizeof(int)) )) return -1;
	for (c=1; c<C; c++) combs[c] = combs[c-1] + d;
	
	/* Define n sets according to the day of the year (julian's day in seismology). */
	for (itr=0; itr<mtr; itr++) {
		tm = gmtime (time + itr);
		setn[itr] = (unsigned int)floor((double)(tm->tm_yday * n)/365.);   /* (0 - 365)/365 */
	}
	
	/* List all combinations */
	for (i=0; i<d; i++) combs[0][i] = i;
	for (c=1; c<C; c++) {
		for (i=0; i<d; i++) combs[c][i] = combs[c-1][i];
		if (combs[c][d-1] < n-1) combs[c][d-1] += 1; 
		else {
			for (i=2; i<=d; i++) 
				if (combs[c][d-i] < n-i) break;
			combs[c][d-i] += 1;
			for (j=d-i+1; j<d; j++) combs[c][j] = combs[c][j-1] + 1;
		}
	}
	
	/* Mark whether a traces is used or not on each combination. */
	for (c=0; c<C; c++) 
		for (itr=0; itr<mtr; itr++) {
			ua1 = setn[itr];
			selC[c][itr] = 1;
			for (i=0; i<d; i++) 
				if (ua1 == combs[c][i]) selC[c][itr] = 0;
		}
	
	myfree(combs[0]);
	myfree(combs);
	myfree(setn);
	
	return 0;
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
int tspws_stacks_float (t_CWTvar *pWTst, t_CWTvar *pWTps, float *sigall, double *w, int max, 
			unsigned int nsamp, unsigned int mtr, t_WaveletFamily *pWF, t_CWTvar *pWT) {
	double complex ca1;
	double da1, *dsig;
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
			
			da1 = (w != NULL) ? w[itr] : 1;
			#pragma omp parallel for default(shared) private(n, ca1) schedule(dynamic,1)
			for (s=0; s<pWT->S; s++) {
				for (n=0; n<pWT->N[s]; n++) {
					ca1 = pWT->d[s][n];
					pWTst->d[s][n] += da1*ca1;            /* Linear stack. */
					ca1 /= cabs(ca1);
					if (creal(ca1*conj(ca1)) <= 1.001) pWTps->d[s][n] += ca1;   /* Phase stack.  */
				}
			}
		}
		myfree(dsig);
	}
	return 0;
}

int tspws_subsmpl_float (float **ls_out, float **tsPWS_out, float *sigall, double *w, int max, unsigned int nsamp, 
			unsigned int mtr, unsigned int K, unsigned int M, double wu, int unbiased, t_WaveletFamily *pWF) {
	t_CWTvar *pWT, *pWTst, *pWTps, *pWTout;
	double complex ca1;
	double *dsig, *W;
	unsigned int m, itr, s, n, nerr = 0;
	float fa1, *sig, *pf1;
	double da1;
	char *sel[M];
	
	if (NULL == (pWT = CreateComplexWaveletVar (pWF, nsamp) )) nerr = 4; /* The trace in time-scale         */
	if (NULL == (pWTst = CreateComplexWaveletVarArray (pWF, nsamp, M) )) nerr = 4; /* The stacked part (& the result) */
	if (NULL == (pWTps = CreateComplexWaveletVarArray (pWF, nsamp, M) )) nerr = 4; /* The phase weight part.          */
	if (NULL == (pWTout = CreateComplexWaveletVar (pWF, nsamp) )) nerr = 4; /* The phase weight part.          */
	if (NULL == (dsig = (double *)mymalloc(max * sizeof(double)) )) nerr = 4;
	if (NULL == (sel[0] = (char *)mymalloc(M*mtr*sizeof(char)) )) nerr = 4;
	if (NULL == (W = (double *)mymalloc(M * sizeof(double)) )) nerr = 4;
	
	if (!nerr) {
		/* Subsampling initialization. */
		for (m=1; m<M; m++) sel[m] = sel[m-1] + mtr;
		for (m=0; m<M; m++) SubsamplingPlan (sel[m], mtr, K);
		for (m=0; m<M; m++) {
			memset(ls_out[m],    0, max*sizeof(float));
			memset(tsPWS_out[m], 0, max*sizeof(float));
		}
		
		for (itr=0; itr<mtr; itr++) {
			sig = sigall + itr*max;
			
			/* WT */
			for(n=0; n<nsamp; n++) dsig[n] = (double)sig[n];
			for(   ; n<max;   n++) dsig[n] = 0.;
			complex_1D_wavelet_dec (pWT, dsig, nsamp, pWF);
			
			/* Linear stacks */
			da1 = (w != NULL) ? w[itr] : 1;
			for (m=0; m<M; m++) {
				pf1 = ls_out[m];
				if (sel[m][itr] == 1)
					for (n=0; n<nsamp; n++) pf1[n] += da1*dsig[n];
			}
			
			/* ts-PWS */ 
			#pragma omp parallel
			{
				/* Linear stack */
				#pragma omp for private(s, n) schedule(dynamic,1)
				for (m=0; m<M; m++)
					if (sel[m][itr] == 1)
						for (s=0; s<pWT->S; s++)
							for (n=0; n<pWT->N[s]; n++)
								pWTst[m].d[s][n] += da1*pWT->d[s][n];
				/* Phase stack.  */
				#pragma omp for private(n, ca1, m) schedule(dynamic,1)
				for (s=0; s<pWT->S; s++)
					for (n=0; n<pWT->N[s]; n++) {
						ca1  = pWT->d[s][n];
						ca1 /= cabs(ca1);
						if (creal(ca1*conj(ca1)) <= 1.001) 
							for (m=0; m<M; m++)
								if (sel[m][itr] == 1) pWTps[m].d[s][n] += ca1;
					}
			}
		}
		
		/* Out of the loop part of the linear & phase-weighed stack. */
		/* Weight corrections. */
		if (w != NULL) {
			for (m=0; m<M; m++) {
				da1 = 0;
				for (itr=0; itr<mtr; itr++)
					if (sel[m][itr] == 1) da1 += w[itr];
				W[m] = (double)K/da1;
			}
		} else for (m=0; m<M; m++) W[m] = 1.;
		
		/* Linear stack. */
		for (m=0; m<M; m++) {
			fa1 = W[m]/K;
			sig = ls_out[m];
			for (n=0; n<max; n++) sig[n] *= fa1;
		}
		
		/* The actual ts-PWS. */
		for (m=0; m<M; m++) {
			if (wu == 2 && unbiased) 
				tspws_unbiased (pWTout, &pWTst[m], &pWTps[m], nsamp, K, K);
			else tspws_biased (pWTout, &pWTst[m], &pWTps[m], nsamp, K, K, wu);
			/* ICWT */
			Re_complex_1D_wavelet_rec (dsig, pWTout, nsamp, pWF);
			/* Save the results */
			sig = tsPWS_out[m];
			for (n=0; n<nsamp; n++) sig[n] = (float)(W[m]*dsig[n]);
			for (   ; n<max;   n++) sig[n] = 0.;
		}
	} else 
		printf ("Error %d in tspws_subsmpl_float()\n", nerr);
	
	
	myfree(W);
	myfree(sel[0]);
	myfree(dsig);
	DestroyComplexWaveletVar (pWTout);
	DestroyComplexWaveletVarArray (pWTps, M);
	DestroyComplexWaveletVarArray (pWTst, M);
	DestroyComplexWaveletVar (pWT);
	
	return nerr;
}

int TwoStage_subsmpl_float (float **ls_out, float **tsPWS_out, float *sigall, double *w, int max, unsigned int nsamp, 
			unsigned int mtr, unsigned int Kmax, unsigned int K, unsigned int M, double wu, int unbiased, 
			t_WaveletFamily *pWF) {
	int nerr = 0;
	
	#pragma omp parallel reduction(+: nerr)
	{
		t_CWTvar *pWT, *pWTst, *pWTps, *pWTout;
		double *dsig, *dsigN, da1, *pd1, W;
		unsigned int m, itr, gtr, k, n;
		float *sig;
		char *sel;
	
		#pragma omp critical
		{
			if (NULL == (pWT = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1; /* The trace in time-scale         */
			if (NULL == (pWTst = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1; /* The stacked part (& the result) */
			if (NULL == (pWTps = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1; /* The phase weight part.          */
			if (NULL == (pWTout = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1 ; /* The phase weight part.          */
			if (NULL == (dsig = (double *)mymalloc(max*sizeof(double)) )) nerr += 1;
			if (NULL == (dsigN = (double *)mymalloc(Kmax*max*sizeof(double)) )) nerr += 1;
			if (NULL == (sel = (char *)mymalloc(mtr*sizeof(char)) )) nerr += 1;
		}
		
		if (!nerr) {
			#pragma omp for schedule(static) 
			for (m=0; m<M; m++) {
				/* Subsampling initialization. */
				SubsamplingPlan (sel, mtr, K);
				
				/* Partial linear_stacks */
				memset(dsigN, 0, Kmax*max*sizeof(double));
				k = 0;
				for (itr=0; itr<mtr; itr++) {
					sig = sigall + itr*max;
					if (sel[itr] == 1) {
						gtr = (unsigned int)floor((double)(k*Kmax)/(double)K);
						pd1 = dsigN + gtr * max;
						da1 = (w != NULL) ? w[itr] : 1;
						for (n=0; n<nsamp; n++) pd1[n] += da1 * (double)sig[n];
						k++;
					}
				}
				
				/* Linear and phase stack of the partial stacks. */
				tspws_stacks_double (pWTst, pWTps, dsigN, max, nsamp, Kmax, pWF, pWT);
				
				/* Out of the loop part of the phase-weighed stack. */
				if (wu == 2 && unbiased) 
					tspws_unbiased (pWTout, pWTst, pWTps, nsamp, Kmax, K);
				else tspws_biased (pWTout, pWTst, pWTps, nsamp, Kmax, K, wu);
				/* ICWT */
				Re_complex_1D_wavelet_rec (dsig, pWTout, nsamp, pWF);
				
				/* Weight corrections. */
				if (w != NULL) {
					da1 = 0;
					for (itr=0; itr<mtr; itr++)
						if (sel[itr] == 1) da1 += w[itr];
					W = (double)K/da1;
				} W = 1.;
			
				/* Save the results */
				sig = tsPWS_out[m];
				for (n=0; n<nsamp; n++) sig[n] = (float)(W*dsig[n]);
				for (   ; n<max;   n++) sig[n] = 0.;
				
				/* Linear stack */
				memcpy(dsig, dsigN, max*sizeof(double));
				for (k=1; k<Kmax; k++) {
					pd1 = dsigN + k*max;
					for (n=0; n<nsamp; n++) dsig[n] += pd1[n];
				}
				/* Normalization */
				da1 = W/K;
				for (n=0; n<nsamp; n++) dsig[n] *= da1;
				/* Save the results */
				sig = ls_out[m];
				for (n=0; n<nsamp; n++) sig[n] = (float)dsig[n];
				for (   ; n<max;   n++) sig[n] = 0.;
			}
		} else 
			printf ("Error %d in TwoStage_subsmpl_float()\n", nerr);
		
		#pragma omp critical
		{
			myfree(sel);
			myfree(dsigN);
			myfree(dsig);
			DestroyComplexWaveletVar (pWTout);
			DestroyComplexWaveletVar (pWTps);
			DestroyComplexWaveletVar (pWTst);
			DestroyComplexWaveletVar (pWT);
		}
	}
	
	return nerr;
}

int tspws_jackknife_float (float **ls_out, float **tsPWS_out, unsigned int *mtr_out, float *sigall, time_t *time, double *w, 
	int max, unsigned int nsamp, unsigned int mtr, unsigned int K, unsigned int d, unsigned int M, 
	unsigned int C, double wu, int unbiased, t_WaveletFamily *pWF) {
	
	return 0;
}

/* Low-memory version of jackknife, this is the similar to the subsmpl version but changing SubsamplingPlan(). */
int TwoStage_jackknife_float (float **ls_out, float **tsPWS_out, unsigned int *mtr_out, float *sigall, time_t *time, double *w, 
	int max, unsigned int nsamp, unsigned int mtr, unsigned int Kmax, unsigned int K0, unsigned int d, 
	unsigned int M, unsigned int C, double wu, int unbiased, t_WaveletFamily *pWF) {
	
	int nerr = 0;
	unsigned int c;
	char **selC;
	
	if (NULL == (selC = (char **)mymalloc(C*sizeof(char *)) )) nerr += 1;
	if (NULL == (selC[0] = (char *)mymalloc(C*mtr*sizeof(char)) )) nerr += 1;
	if (nerr) return nerr;
	for (c=1; c<C; c++) selC[c] = selC[c-1] + mtr;
	
	/* Subsampling initialization. */
	JackknifePlans (selC, time, mtr, d, M, C);
	
	#pragma omp parallel default(shared) reduction(+: nerr)
	{
		t_CWTvar *pWT, *pWTst, *pWTps, *pWTout;
		double *dsig, *dsigN, da1, *pd1, W;
		unsigned int m, itr, gtr, k, K, n;
		float *sig;
		char *sel;
		
		#pragma omp critical
		{
			if (NULL == (pWT = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1; /* The trace in time-scale         */
			if (NULL == (pWTst = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1; /* The stacked part (& the result) */
			if (NULL == (pWTps = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1; /* The phase weight part.          */
			if (NULL == (pWTout = CreateComplexWaveletVar (pWF, nsamp) )) nerr += 1 ; /* The phase weight part.          */
			if (NULL == (dsig = (double *)mymalloc(max*sizeof(double)) )) nerr += 1;
			if (NULL == (dsigN = (double *)mymalloc(Kmax*max*sizeof(double)) )) nerr += 1;
		}
		
		if (!nerr) {
			#pragma omp for schedule(static) 
			for (m=0; m<C; m++) {
				sel = selC[m];
				
				/* Partial linear_stacks */
				memset(dsigN, 0, Kmax*max*sizeof(double));
				K = 0;
				for (itr=0; itr<mtr; itr++) if (sel[itr] == 1) K++;
				k = 0;
				for (itr=0; itr<mtr; itr++) {
					sig = sigall + itr*max;
					if (sel[itr] == 1) {
						gtr = (unsigned int)floor((double)(k*Kmax)/(double)K);
						pd1 = dsigN + gtr * max;
						da1 = (w != NULL) ? w[itr] : 1;
						for (n=0; n<nsamp; n++) pd1[n] += da1 * (double)sig[n];
						k++;
					}
				}
				
				/* Linear and phase stack of the partial stacks. */
				tspws_stacks_double (pWTst, pWTps, dsigN, max, nsamp, Kmax, pWF, pWT);
				
				/* Out of the loop part of the phase-weighed stack. */
				if (wu == 2 && unbiased) 
					tspws_unbiased (pWTout, pWTst, pWTps, nsamp, Kmax, K);
				else tspws_biased (pWTout, pWTst, pWTps, nsamp, Kmax, K, wu);
				/* ICWT */
				Re_complex_1D_wavelet_rec (dsig, pWTout, nsamp, pWF);
				
				/* Weight corrections. */
				if (w != NULL) {
					da1 = 0;
					for (itr=0; itr<mtr; itr++)
						if (sel[itr] == 1) da1 += w[itr];
					W = (double)K/da1;
				} W = 1.;
			
				/* Save the results */
				mtr_out[m] = K;
				
				sig = tsPWS_out[m];
				for (n=0; n<nsamp; n++) sig[n] = (float)(W*dsig[n]);
				for (   ; n<max;   n++) sig[n] = 0.;
				
				/* Linear stack */
				memcpy(dsig, dsigN, max*sizeof(double));
				for (k=1; k<Kmax; k++) {
					pd1 = dsigN + k*max;
					for (n=0; n<nsamp; n++) dsig[n] += pd1[n];
				}
				/* Normalization */
				da1 = W/K;
				for (n=0; n<nsamp; n++) dsig[n] *= da1;
				/* Save the results */
				sig = ls_out[m];
				for (n=0; n<nsamp; n++) sig[n] = (float)dsig[n];
				for (   ; n<max;   n++) sig[n] = 0.;
			}
		} else 
			printf ("Error %d in TwoStage_subsmpl_float()\n", nerr);
		
		#pragma omp critical
		{
			myfree(dsigN);
			myfree(dsig);
			DestroyComplexWaveletVar (pWTout);
			DestroyComplexWaveletVar (pWTps);
			DestroyComplexWaveletVar (pWTst);
			DestroyComplexWaveletVar (pWT);
		}
	}
	
	myfree(selC[0]);
	myfree(selC);
	
	return nerr;
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
void partial_linear_stacks (double *dsigN, float *sigall, double *w, unsigned int max, unsigned int nsamp, 
			unsigned int Kmax, unsigned int mtr) {
	unsigned int itr, gtr, n;
	float *sig;
	double da1, *pd1;
	
	memset(dsigN, 0, Kmax*max*sizeof(double));
	for (itr=0; itr<mtr; itr++) {
		sig = sigall + itr*max;
		/* gtr = (itr % Kmax) */  /* Interlaced version */
		gtr = (unsigned int)floor((double)(itr*Kmax)/(double)mtr);
		pd1 = dsigN + gtr*max;
		da1 = (w != NULL) ? w[itr] : 1;
		for (n=0; n<nsamp; n++) pd1[n] += da1 * (double)sig[n];
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
		
		// #pragma omp parallel for default(shared) private(n, ca1) schedule(dynamic,1)
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
		da1 = 1./(da1*da1*(double)mtr);
		// #pragma omp parallel for default(shared) private(n, ca1, da2) schedule(dynamic,1)
		for (s=0; s<pWTst->S; s++) {
			for (n=0; n<pWTst->N[s]; n++) {
				ca1 = pWTps->d[s][n];
				da2 = (creal(ca1) * creal(ca1) + cimag(ca1) * cimag(ca1)) * da1;
				pWTout->d[s][n] = da2 * pWTst->d[s][n];
			}
		}
	} else if (wu == 1) {
		da1 = 1./((double)Kmax*(double)mtr);
		// #pragma omp parallel for default(shared) private(n) schedule(dynamic,1)
		for (s=0; s<pWTst->S; s++)
			for (n=0; n<pWTst->N[s]; n++)
				pWTout->d[s][n] = pWTst->d[s][n] * cabs(pWTps->d[s][n]) * da1;
	} else {
		// #pragma omp parallel for default(shared) private(n, ca1, da2) schedule(dynamic,1)
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
#if 0
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
#else
void tspws_unbiased (t_CWTvar *pWTout, t_CWTvar *pWTst, t_CWTvar *pWTps, unsigned int nsamp, 
			unsigned int Kmax, unsigned int mtr) {
	double complex ca1;
	double da1;
	double invKmax = 1./(double)Kmax, invKmax_1 = 1./(double)(Kmax-1), invmtr = 1./(double)mtr;
	unsigned int n, s;

	if (Kmax == 1) tspws_biased (pWTout, pWTst, pWTps, nsamp, Kmax, mtr, 2);
	else {
		// #pragma omp parallel for default(shared) private(n, ca1, da1) schedule(dynamic,1)
		for (s=0; s<pWTst->S; s++) {
			for (n=0; n<pWTst->N[s]; n++) {
				ca1 = pWTps->d[s][n] * invKmax;
				da1 = creal(ca1) * creal(ca1) + cimag(ca1) * cimag(ca1);
				da1 = (Kmax*da1 - 1) * invKmax_1;
				pWTout->d[s][n] = pWTst->d[s][n] * da1 * invmtr;
			}
		}
	}
}
#endif
