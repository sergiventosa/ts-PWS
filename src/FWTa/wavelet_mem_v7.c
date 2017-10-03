/*
 * Minimal wavelet library using frames. Based on the wavelet library developed in my PhD. for the time-scale slowness filters.
 * - Memory routines.
 * 2011/01 Sergi Ventosa
 * 2016/01 Sergi Ventosa: migrated to c99 double complex types.
 */
#include <string.h>

#include "wavelet_v7.h"
#include "myallocs.h"
#include "prnmsg.h"

/**********************************/
/* Create a real wavelet variable */
/**********************************/
t_RWTvar *Create_RWV_error(t_RWTvar *p) {
	prerror("CreateRealWaveletVar run out of memory.");
	DestroyRealWaveletVar(p);
	return p;
}

t_CWTvar *Create_CWV_error(t_CWTvar *p) {
	prerror("CreateComplexWaveletVar run out of memory.");
	DestroyComplexWaveletVar(p);
	return p;
}

t_RWTvar *Create_RWVA_error(t_RWTvar *p, unsigned int sz) {
	prerror("CreateRealWaveletVarArray run out of memory.");
	DestroyRealWaveletVarArray(p, sz);
	return p;
}

t_CWTvar *Create_CWVA_error(t_CWTvar *p, unsigned int sz) {
	prerror("CreateComplexWaveletVarArray run out of memory.");
	DestroyComplexWaveletVarArray(p, sz);
	return p;
}

int CheckWaveletFamily (t_WaveletFamily *pWF) {
	if (pWF == NULL || pWF->Down_smp == NULL || pWF->Ns <= 0) {
		prerror("The waveletFamily variable is not nice.");
		return 1;
	}
	return 0;
}

t_RWTvar *CreateRealWaveletVar (t_WaveletFamily *pWF, unsigned int N) {
	t_RWTvar *p;
	unsigned int s, S, ne, ia1;

	if (CheckWaveletFamily (pWF)) return NULL;
	S = pWF->Ns;
	if (NULL == (p = (t_RWTvar *)mycalloc(1, sizeof(t_RWTvar)) ))
		return Create_RWV_error(p);
	if (NULL == (p->N = (unsigned int *)mycalloc(S, sizeof(int)) ))
		return Create_RWV_error(p);
	if (NULL == (p->d = (double **)mycalloc(S, sizeof(double *)) ))
		return Create_RWV_error(p);
	p->S = S;

	ne = 0;
	for (s=0; s<S; s++) {
		ia1 = (N + pWF->Down_smp[s] -1)/pWF->Down_smp[s];
		ne += ia1;
		p->N[s] = ia1;
	}

	if (NULL == (p->d[0] = (double *)mycalloc(ne, sizeof(double)) ))
		return Create_RWV_error(p);
		
	for (s=0; s<S-1; s++) 
		p->d[s+1] = p->d[s] + p->N[s];

	return p;
}

t_CWTvar *CreateComplexWaveletVar (t_WaveletFamily *pWF, unsigned int N) {
	t_CWTvar *p;
	unsigned int s, S, ia1, ia2;

	if (CheckWaveletFamily (pWF)) return NULL;
	S = pWF->Ns;
	if (NULL == (p = (t_CWTvar *)mycalloc(1, sizeof(t_CWTvar)) ))
		return Create_CWV_error(p);
	if (NULL == (p->N = (unsigned int *)mycalloc(S, sizeof(int)) ))
		return Create_CWV_error(p);
	if (NULL == (p->d = (double complex **)mycalloc(S, sizeof(double complex *)) ))
		return Create_CWV_error(p);
	p->S = S;

	ia2 = 0;
	for (s=0; s<S; s++) {
		ia1 = (N + pWF->Down_smp[s] -1)/pWF->Down_smp[s];
		ia2 += ia1;
		p->N[s] = ia1;
	}

	if (NULL == (p->d[0] = (double complex *)mycalloc(ia2, sizeof(double complex)) ))
		return Create_CWV_error(p);
		
	for (s=0; s<S-1; s++) 
		p->d[s+1] = p->d[s] + p->N[s];

	return p;
}

t_RWTvar *CreateRealWaveletVarArray (t_WaveletFamily *pWF, unsigned int N, unsigned int sz)  {
	t_RWTvar *p, *p0;
	unsigned int i, s, S, ne, ia1;

	if (CheckWaveletFamily (pWF)) return NULL;
	S = pWF->Ns;
	if (sz <= 0) {prerror("CreateRealWaveletVarArray"); return NULL; }
	if (NULL == (p = (t_RWTvar *)mycalloc(sz, sizeof(t_RWTvar)) ))
		return Create_RWVA_error(p, sz);
	if (NULL == (p->N = (unsigned int *)mycalloc(sz*S, sizeof(int)) ))
		return Create_RWVA_error(p, sz);
	if (NULL == (p->d = (double **)mycalloc(sz*S, sizeof(double *)) ))
		return Create_RWVA_error(p, sz);
	p->S = S;

	ne = 0;
	for (s=0; s<S; s++) {
		ia1 = (N + pWF->Down_smp[s] -1)/pWF->Down_smp[s];
		ne += ia1;
		p->N[s] = ia1;
	}

	if (NULL == (p->d[0] = (double *)mycalloc(sz*ne, sizeof(double)) ))
		return Create_RWVA_error(p, sz);
		
	for (s=0; s<S-1; s++) 
		p->d[s+1] = p->d[s] + p->N[s];

	p0 = p;
	for (i=1; i<sz; i++) {
		p++;
		p->S = S;
		p->N = p0->N + S*i;
		p->d = p0->d + S*i;
		memcpy (p->N, p0->N, S*sizeof(int));
		for (s=0; s<S; s++)
			p->d[s] = p0->d[s] + ne*i;
	}

	return p0;
}

t_CWTvar *CreateComplexWaveletVarArray (t_WaveletFamily *pWF, unsigned int N, unsigned int sz)  {
	t_CWTvar *p, *p0;
	unsigned int i, s, S, ne, ia1;

	if (CheckWaveletFamily (pWF)) return NULL;
	S = pWF->Ns;
	if (sz <= 0) {prerror("CreateComplexWaveletVarArray"); return NULL; }
	if (NULL == (p = (t_CWTvar *)mycalloc(sz, sizeof(t_CWTvar)) ))
		return Create_CWVA_error(p, sz);
	if (NULL == (p->N = (unsigned int *)mycalloc(sz*S, sizeof(int)) ))
		return Create_CWVA_error(p, sz);
	if (NULL == (p->d = (double complex **)mycalloc(sz*S, sizeof(double complex *)) ))
		return Create_CWVA_error(p, sz);
	p->S = S;

	ne = 0;
	for (s=0; s<S; s++) {
		ia1 = (N + pWF->Down_smp[s] -1)/pWF->Down_smp[s];
		ne += ia1;
		p->N[s] = ia1;
	}

	if (NULL == (p->d[0] = (double complex *)mycalloc(sz*ne, sizeof(double complex)) ))
		return Create_CWVA_error(p, sz);
		
	for (s=0; s<S-1; s++) 
		p->d[s+1] = p->d[s] + p->N[s];

	p0 = p;
	for (i=1; i<sz; i++) {
		p++;
		p->S = S;
		p->N = p0->N + S*i;
		p->d = p0->d + S*i;
		memcpy (p->N, p0->N, S*sizeof(int));
		for (s=0; s<S; s++)
			p->d[s] = p0->d[s] + ne*i;
	}

	return p0;
}

void CleanRealWaveletVar (t_RWTvar *p) {
	unsigned int s, S, ne;
	
	if (p != NULL) {
		S = p->S;
		ne = 0;
		for (s=0; s<S; s++) ne += p->N[s];
		memset(p->d[0], 0, ne*sizeof(double));
	}
}

void CleanComplexWaveletVar (t_CWTvar *p) {
	unsigned int s, S, ne;
	
	if (p != NULL) {
		S = p->S;
		ne = 0;
		for (s=0; s<S; s++) ne += p->N[s];
		memset(p->d[0], 0, ne*sizeof(double complex));
	}
}

void CleanRealWaveletVarArray (t_RWTvar *p, unsigned int sz) {
	unsigned int s, S, ne;
	
	if (p != NULL) {
		S = p->S;
		ne = 0;
		for (s=0; s<S; s++) ne += p->N[s];
		memset(p->d[0], 0, sz*ne*sizeof(double));
	}
}

void CleanComplexWaveletVarArray (t_CWTvar *p, unsigned int sz) {
	unsigned int s, S, ne;
	
	if (p != NULL) {
		S = p->S;
		ne = 0;
		for (s=0; s<S; s++) ne += p->N[s];
		memset(p->d[0], 0, sz*ne*sizeof(double complex));
	}
}

void DestroyRealWaveletVar (t_RWTvar *p) {
	if (p != NULL) {
		myfree(p->d[0]);
		myfree(p->d);
		myfree(p->N);
		myfree(p);
		p = NULL;
	}
}

void DestroyComplexWaveletVar (t_CWTvar *p) {
	if (p != NULL) {
		myfree(p->d[0]);
		myfree(p->d);
		myfree(p->N);
		myfree(p);
		p = NULL;
	}
}

void DestroyRealWaveletVarArray (t_RWTvar *p, unsigned int sz) {
	DestroyRealWaveletVar (p);
}

void DestroyComplexWaveletVarArray (t_CWTvar *p, unsigned int sz) {
	DestroyComplexWaveletVar (p);
}

