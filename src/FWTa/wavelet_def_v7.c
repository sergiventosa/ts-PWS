/*
 * Minimal wavelet library using frames. Based on the wavelet library developed in my PhD. for the time-scale slowness filters.
 * - Wavelet function definitions.
 * 2011/01 Sergi Ventosa
 * 2016/01 Sergi Ventosa: migrated to c99 double complex types.
 * 2018/04 Sergi Ventosa: Change normalization of MexHat, Cpsi & IFWT to be close to the most commontly definitions.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "wavelet_v7.h"
#include "myallocs.h"
#include "prnmsg.h"

void setscales0 (double *ps, double s0, unsigned int S, unsigned int V, double a0);
unsigned int setwaveletlength0 (unsigned int *Ls, double *ps, unsigned int S, unsigned int MaxL);
void setsampling0 (unsigned int *DownSmpl, unsigned int S, int continuous, unsigned int V, unsigned int J, double b0, double a0, double s0);

/* Build a family of Morlet wavelets at the specified scales. */
/* (Approximated version, the non-zero mean one) */
void MorletFun_real (double *w, int *center, double w0, unsigned int Ls, double scale) {
	double k, daux, w1;
	int u, u0;

	k = 1/sqrt(sqrt(PI)*scale);
	*center = Ls/2;
	scale = 1/scale;
	u0 = Ls/2;
	for (u = -u0; u < (signed)Ls - u0; u++) {
		daux = scale * u;
		w1 = w0*daux;
		daux *= daux;
		*w++ = k*cos(w1)*exp(-0.5*daux);
	}
}

void MorletFun (double complex *w, int *center, double w0, unsigned int Ls, double scale) {	
	double k, daux, w1;
	int u, u0;

	k = 1/sqrt(sqrt(PI)*scale);
	*center = Ls/2;
	scale = 1/scale;
	u0 = Ls/2;
	for (u = -u0; u < (signed)Ls - u0; u++) {
		daux = scale * u;
		w1 = w0*daux;
		daux *= daux;
		*w++ = k * cexp(I*w1) * exp(-0.5*daux);
	}
}

void Complete_MorletFun_real (double *w, int *center, double w0, unsigned int Ls, double scale) {
	const double ze = exp((-w0*w0)/2);
	double k, daux, w1;
	int u, u0;

	k = 1/sqrt(sqrt(PI)*scale);
	*center = Ls/2;
	scale = 1/scale;
	u0 = Ls/2;
	for (u = -u0; u < (signed)Ls - u0; u++) {
		daux = scale * u;
		w1 = w0*daux;
		daux *= daux;
		*w++ = k*(cos(w1)-ze)*exp(-0.5*daux);
	}
}

void Complete_MorletFun (double complex *w, int *center, double w0, unsigned int Ls, double scale) {
	const double ze = exp((-w0*w0)/2);
	double k, daux, w1;
	int u, u0;

	k = 1/sqrt(sqrt(PI)*scale);
	*center = Ls/2;
	scale = 1/scale;
	u0 = Ls/2;
	for (u = -u0; u < (signed)Ls - u0; u++) {
		daux = scale * u;
		w1 = w0*daux;
		daux *= daux;
		daux = k * exp(-0.5*daux);
		*w++ = daux * (cexp(I*w1) - ze);
	}
}

void MexicanHatFun_real (double *w, int *center, double dump, unsigned int Ls, double scale) {
	double k, daux;
	int u, u0;

	k = 2/sqrt(3*sqrt(PI)*scale);
	*center = Ls/2;
	scale = 1/scale;
	u0 = Ls/2;
	for (u = -u0; u < (signed)Ls - u0; u++) {
		daux  = scale*u;
		daux *= daux;
		*w++  = k *(daux - 1)*exp(-daux/2);
	}
}

double erfi (double z) {
	double da1, out, zz = z*z;
	unsigned int n;
	
	da1 = z;
	out = da1;
	for (n=1; n<500; n++) {
		da1 *= zz / n;
		out += da1 / (2*n+1);
	}
	out *= 2/sqrt(PI);
	
	return out;
}

void MexicanHatFun (double complex *w, int *center, double dump, unsigned int Ls, double scale) {
	double k, daux;
	int u, u0;

	k = 2/sqrt(3*sqrt(PI)*scale);
	*center = Ls/2;
	scale = 1/scale;
	u0 = Ls/2;
	for (u = -u0; u < (signed)Ls - u0; u++) {
		daux = scale*u;
		*w++ = k * ( (daux*daux - 1)*exp(-daux*daux/2) * (1 + I*erfi(daux/sqrt(2))) - I*sqrt(2/PI)*daux);
	}
}

double Morlet_Cpsi (double w0) {
	double w, Cpsi=0, aux;

	for (w=0.01; w<100; w+=0.01) {
		aux = (w-w0);
		aux *= aux;
		Cpsi += exp(-aux)/w;
	}
	Cpsi *= 0.01*sqrt(PI)/2;

	return (Cpsi);
}	

double MexicanHat_Cpsi () {
	return ( (4. / 3.) * sqrt(PI));
}


/* Same as above. */
int FillDualFrame(t_WaveletFamily *pWF) {
	unsigned int s, n, Nsmpl;
	
	Nsmpl = 0;
	for (s=0; s<pWF->Ns; s++) {
		pWF->Lds[s] = pWF->Ls[s];
		Nsmpl += pWF->Ls[s];
	}
	
	for (s=0; s<pWF->Ns; s++) pWF->center_df[s] = pWF->Ls[s] - 1 - pWF->center[s];

	if (pWF->format == 1) {
		double *din, *dout;
		
		if (pWF->wdualframe.wr[0] != NULL) myfree(pWF->wdualframe.wr[0]);
		if (NULL == (pWF->wdualframe.wr[0] = (double *)mymalloc(Nsmpl*sizeof(double))) ) return 1;
		for (s=1; s<pWF->Ns; s++) pWF->wdualframe.wr[s] = pWF->wdualframe.wr[s-1] + pWF->Lds[s-1];
		for (s=0; s<pWF->Ns; s++) {
			din = pWF->wframe.wr[s] + pWF->Ls[s] - 1;
			dout = pWF->wdualframe.wr[s];
			for (n=0; n<pWF->Ls[s]; n++) *dout++ = *din--;
		}
	} else if (pWF->format == -1) {
		double complex *cin, *cout;
		
		if (pWF->wdualframe.wc[0] != NULL) myfree(pWF->wdualframe.wc[0]);
		if (NULL == (pWF->wdualframe.wc[0] = (double complex *)mymalloc(Nsmpl*sizeof(double complex))) ) return 1;
		for (s=1; s<pWF->Ns; s++) pWF->wdualframe.wc[s] = pWF->wdualframe.wc[s-1] + pWF->Lds[s-1];
		for (s=0; s<pWF->Ns; s++) {
			cin = pWF->wframe.wc[s] + pWF->Ls[s] - 1;
			cout = pWF->wdualframe.wc[s];
			for (n=0; n<pWF->Ls[s]; n++)
				*cout++ = conj(*cin--);
		}
	}
	return 0;
}

t_WaveletFamily *Create_WF_error(t_WaveletFamily *p) {
	prerror("CreateWaveletFamily run out of memory.");
	DestroyWaveletFamily(p);
	return p;
}

/*************************************************************/
/* Create a Morlet wavelet family from a set of Ns scales.   */
/* Inputs:                                                   */
/*  scales: Double vector with the scales.                   */
/*  Ns: Length of the scale vector.                          */
/* Outputs:                                                  */
/*  Pointer to the t_Wavelet type structure created          */
/*************************************************************/
t_WaveletFamily *CreateWaveletFamily (int type, unsigned int Ns, unsigned int V, unsigned int N, double s0, double b0, int convtype, double op1, int continuous) {
	t_WaveletFamily *pWF;
	unsigned int s, Nsmpl;
	
	void (*RealWavelets[])(double *, int *, double, unsigned int, double) = {NULL, MorletFun_real, Complete_MorletFun_real, MexicanHatFun_real};
	void (*ComplexWavelets[])(double complex *, int *, double, unsigned int, double) = {NULL, MorletFun, Complete_MorletFun, MexicanHatFun};

	/* Reserve memory */
	if (!type) {
		prerror("CreateWaveletFamily: Wavelet type is 0.");
		return NULL;
	}
	if (type > 0 && type > MAXREAL_WAVE) {
		prerror("CreateWaveletFamily: Real wavelet not defined.");
		return NULL;
	}
	if (type < 0 && type < -MAXCMPLX_WAVE) {
		prerror("CreateWaveletFamily: Complex wavelet not defined.");
		return NULL;
	}

	Ns *= V;
	pWF = (t_WaveletFamily *)mycalloc(1, sizeof(t_WaveletFamily));
	pWF->scale = (double *)mymalloc(Ns*sizeof(double));
	pWF->Ls = (unsigned int *)mymalloc(Ns*sizeof(int));
	pWF->Lds = (unsigned int *)mycalloc(Ns,sizeof(int));
	pWF->center = (int *)mymalloc(Ns*sizeof(int));
	pWF->center_df = (int *)mymalloc(Ns*sizeof(int));
	pWF->Down_smp = (unsigned int *)mymalloc(Ns*sizeof(int));
	if (pWF == NULL || pWF->scale == NULL || pWF->Ls == NULL || pWF->Lds == NULL || pWF->center == NULL || pWF->center_df == NULL || pWF->Down_smp == NULL) 
		return Create_WF_error(pWF);

	pWF->type = type;
	if (type > 0) {
		pWF->format = 1;
		if (NULL == (pWF->wframe.wr = (double **)mycalloc(Ns,sizeof(double *)) ))
			return Create_WF_error(pWF);
		if (NULL == (pWF->wdualframe.wr = (double **)mycalloc(Ns,sizeof(double *)) ))
			return Create_WF_error(pWF);
	} else if (type <  0) {
		pWF->format = -1;
		if (NULL == (pWF->wframe.wc = (double complex **)mycalloc(Ns,sizeof(double complex *)) ))
			return Create_WF_error(pWF);
		if (NULL == (pWF->wdualframe.wc = (double complex **)mycalloc(Ns,sizeof(double complex *)) ))
			return Create_WF_error(pWF);
	} else return Create_WF_error(pWF);

	/* Fill it */
	pWF->a0 = 2;
	pWF->b0 = continuous ? 1 : b0;
	pWF->convtype = convtype;
	pWF->op1 = op1;
	pWF->V = V;
	pWF->Ns = Ns;
	setscales0 (pWF->scale, s0, Ns, V, 2);
	Nsmpl = setwaveletlength0 (pWF->Ls, pWF->scale, Ns, N);
	setsampling0 (pWF->Down_smp, Ns, continuous, V, Ns/V, pWF->b0, 2, s0);

	if (type > 0) {
		pWF->wdualframe.wr[0] = NULL;
		if (NULL == (pWF->wframe.wr[0] = (double *)mymalloc(Nsmpl*sizeof(double))) )
			return Create_WF_error(pWF);
		for (s=1; s<pWF->Ns; s++) pWF->wframe.wr[s] = pWF->wframe.wr[s-1] + pWF->Ls[s-1];
		for (s=0; s<Ns; s++) {
			(*RealWavelets[type])(pWF->wframe.wr[s], pWF->center + s, op1, pWF->Ls[s], pWF->scale[s]);
		}
	} else if (type < 0) {
		pWF->wdualframe.wc[0] = NULL;
		if (NULL == (pWF->wframe.wc[0] = (double complex *)mymalloc(Nsmpl*sizeof(double complex))) ) 
			return Create_WF_error(pWF);
		for (s=1; s<pWF->Ns; s++) pWF->wframe.wc[s] = pWF->wframe.wc[s-1] + pWF->Ls[s-1];
		for (s=0; s<Ns; s++) 
			(*ComplexWavelets[-type])(pWF->wframe.wc[s], pWF->center + s, op1, pWF->Ls[s], pWF->scale[s]);
	} else return Create_WF_error(pWF);
	
	if (FillDualFrame(pWF)) return Create_WF_error(pWF);

	/* Cpsi */
	switch (type) {
	case -2:
	case -1: 
	case  1:
	case  2: pWF->Cpsi = Morlet_Cpsi (op1);
	         break;
	case -3:
	case  3: pWF->Cpsi = MexicanHat_Cpsi ();
	         break;
	default: pWF->Cpsi = 1;
	}
	
	return pWF;
}


/* Set the values of the scales sampled. */
void setscales0 (double *ps, double s0, unsigned int S, unsigned int V, double a0) {
	double da1, da2;
	unsigned int s;
	
	da1 = s0;
	da2 = pow(a0, 1/(double)V);
	for (s=0; s<S; s++) {
		ps[s] = da1;
		da1 *= da2;
	}
}

/* Set wavelet functions length. */
unsigned int setwaveletlength0 (unsigned int *Ls, double *ps, unsigned int S, unsigned int MaxL) {
	unsigned int s, Nsmpl, ua1;
	
	Nsmpl = 0;
	for (s=0; s<S; s++) {
		ua1 = 2*(unsigned int)ceil(NSIGMAS*ps[s]) + 1;
		Ls[s] = (ua1>MaxL) ? MaxL : ua1;
		Nsmpl += ua1;
	}
	return Nsmpl;
}

/* Set Downsampling field. */
void setsampling0 (unsigned int *DownSmpl, unsigned int S, int continuous, unsigned int V, unsigned int J, double b0, double a0, double s0) {
	double da1;
	unsigned int s, j, v, ua1;
	
	if (continuous) 
		for (s=0; s<S; s++) DownSmpl[s] = 1;
	else {
		da1 = s0*b0;
		for (j=0; j<J; j++) {
			ua1 = (da1 < 1) ? 1 : (unsigned int)da1;
			for (v=0; v<V; v++) *DownSmpl++ = ua1;
			da1 *= a0;
		}
	}
}

/*****************************************************/
/* Free wavelet family memory                        */
/*  pWF: Pointer to t_WaveletFamily type structure.  */
/*****************************************************/
void DestroyWaveletFamily (t_WaveletFamily *pWF) {
	if (pWF != NULL) {
		if (pWF->format == 1) {
			if (pWF->wframe.wr != NULL) myfree(pWF->wframe.wr[0]);
			myfree (pWF->wframe.wr);
			if (pWF->wdualframe.wr != NULL) myfree(pWF->wdualframe.wr[0]);
			myfree (pWF->wdualframe.wr);
		} else if (pWF->format == -1) {
			if (pWF->wframe.wc != NULL) myfree(pWF->wframe.wc[0]);
			myfree (pWF->wframe.wc);
			if (pWF->wdualframe.wc != NULL) myfree(pWF->wdualframe.wc[0]);
			myfree (pWF->wdualframe.wc);
		}
		myfree (pWF->Down_smp);
		myfree (pWF->center_df);
		myfree (pWF->center);
		myfree (pWF->Lds);
		myfree (pWF->Ls);
		myfree (pWF->scale);
		myfree (pWF);
		pWF = NULL;
	}
}
