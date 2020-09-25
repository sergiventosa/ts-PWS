#ifdef MATLAB
#include <string.h>

#include "matrix.h"
#include "MatlabIO.h"
#include "MatlabIO_WT.h"

mxArray *export_RealWaveletDec(double **d, unsigned int *N, unsigned int S) {
	mxArray *pa, *pm;
	unsigned int j;
	mwSize dims[2];

	if (d==NULL || N==NULL) return NULL;
	dims[0] = 1;
	dims[1] = S;
	pa = mxCreateCellArray(2, dims);
	for (j=0; j<S; j++) {
		pm = mxCreateDoubleMatrix(1, N[j], mxREAL);
		mxSetCell(pa, j, pm);
		memcpy(mxGetPr(pm), d[j], N[j]*sizeof(double));
	}
	return pa;
}

mxArray *export_ComplexWaveletDec(double complex **d, unsigned int *N, unsigned int S) {
	mxArray *pa, *pm;
	unsigned int i, j;
	mwSize dims[2];
	double *pr, *pi;
	double complex *wc;

	if (d==NULL || N==NULL) return NULL;
	dims[0] = 1;
	dims[1] = S;
	pa = mxCreateCellArray(2, dims);
	for (j=0; j<S; j++) {
		pm = mxCreateDoubleMatrix(1, N[j], mxCOMPLEX);
		mxSetCell(pa, j, pm);
		pr = mxGetPr(pm);
		pi = mxGetPi(pm);
		wc = d[j];
		for (i=0; i<N[j]; i++) {
			pr[i] = creal(wc[i]);
			pi[i] = cimag(wc[i]);
		}
	}
	return pa;
}

mxArray *export_RWT(t_RWTvar *x, int M, int N) {
	const char *field_names[] = {"d", "N", "S"};
	int m, i, S;
	mwSize dims[2];
	mxArray *pa, *pm;
	double *p;

	dims[0] = M;
	dims[1] = N;
	pa = mxCreateStructArray(2, dims, sizeof(field_names)/sizeof(*field_names), field_names);
	
	if (x==NULL) return NULL;
	for (m=0; m<M*N; m++, x++) {
		S = x->S;
		mxSetField(pa, m, "d", export_RealWaveletDec(x->d, x->N, S));
		
		pm = mxCreateDoubleMatrix(1, S, mxREAL);
		p = mxGetPr(pm);
		for (i=0; i<S; i++) *p++ = (double)x->N[i];
		mxSetField(pa, m, "N", pm);
		
		mxSetField(pa, m, "S", mxCreateDoubleScalar((double)S));
	}
	return pa;
}

mxArray *export_CWT(t_CWTvar *x, int M, int N) {
	const char *field_names[] = {"d", "N", "S"};
	int m, i, S;
	mwSize dims[2];
	mxArray *pa, *pm;
	double *p;

	if (x==NULL) return NULL;
	dims[0] = M;
	dims[1] = N;
	pa = mxCreateStructArray(2, dims, sizeof(field_names)/sizeof(*field_names), field_names);
	
	for (m=0; m<M*N; m++, x++) {
		S = x->S;
		mxSetField(pa, m, "d", export_ComplexWaveletDec(x->d, x->N, S));
		
		pm = mxCreateDoubleMatrix(1, S, mxCOMPLEX);
		p = mxGetPr(pm);
		for (i=0; i<S; i++) *p++ = (double)x->N[i];
		mxSetField(pa, m, "N", pm);
		
		mxSetField(pa, m, "S", mxCreateDoubleScalar((double)S));
	}
	return pa;
}

mxArray *export_WaveletFamily(t_WaveletFamily *pWF) {
	return export_WaveletFamilyArray(&pWF, 1, 1);
}

mxArray *export_WaveletFamilyArray(t_WaveletFamily **pWF, int M, int N) {
	const char *field_names[] = {"format", "type", "convtype", "wave", "scale", "Ls", 
								 "center", "Down_smp", "Ns", "Cpsi", "a0", "b0", "op1", "V"};
	int m, i, Ns;
	mwSize dims[2];
	mxArray *pa, *pm;
	double *p;
	
	if (pWF==NULL) return NULL;
	dims[0] = M;
	dims[1] = N;
	pa = mxCreateStructArray(2, dims, sizeof(field_names)/sizeof(*field_names), field_names);
	
	for (m=0; m<M*N; m++) {
		if (pWF[m] == NULL) continue;
		Ns = pWF[m]->Ns;
		mxSetField(pa, m, "format", mxCreateDoubleScalar((double)pWF[m]->format));
		mxSetField(pa, m, "type", mxCreateDoubleScalar((double)pWF[m]->type));
		mxSetField(pa, m, "convtype", mxCreateDoubleScalar((double)pWF[m]->convtype));

		if (pWF[m]->format < 0) pm = export_ComplexWaveletDec(pWF[m]->wframe.wc, pWF[m]->Ls, Ns);
		else pm = export_RealWaveletDec(pWF[m]->wframe.wr, pWF[m]->Ls, Ns);
		mxSetField(pa, m, "wave", pm);
		
		pm = mxCreateDoubleMatrix(1, Ns, mxREAL);
		memcpy(mxGetPr(pm), pWF[m]->scale, Ns*sizeof(double));
		mxSetField(pa, m, "scale", pm);

		pm = mxCreateDoubleMatrix(1, Ns, mxREAL);
		p = mxGetPr(pm);
		for (i=0; i<Ns; i++) *p++ = (double)pWF[m]->Ls[i];
		mxSetField(pa, m, "Ls", pm);

		pm = mxCreateDoubleMatrix(1, Ns, mxREAL);
		p = mxGetPr(pm);
		for (i=0; i<Ns; i++) *p++ = (double)pWF[m]->center[i];
		mxSetField(pa, m, "center", pm);

		pm = mxCreateDoubleMatrix(1, Ns, mxREAL);
		p = mxGetPr(pm);
		for (i=0; i<Ns; i++) *p++ = (double)pWF[m]->Down_smp[i];
		mxSetField(pa, m, "Down_smp", pm);

		mxSetField(pa, m, "Ns", mxCreateDoubleScalar((double)Ns));
		mxSetField(pa, m, "Cpsi", mxCreateDoubleScalar(pWF[m]->Cpsi));
		mxSetField(pa, m, "a0", mxCreateDoubleScalar(pWF[m]->a0));
		mxSetField(pa, m, "b0", mxCreateDoubleScalar(pWF[m]->b0));
		mxSetField(pa, m, "op1", mxCreateDoubleScalar(pWF[m]->op1));
		mxSetField(pa, m, "V", mxCreateDoubleScalar((double)pWF[m]->V));
	}
	return pa;
}

#endif

