#ifdef MATLAB
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "wavelet_v7.h"
#include "MatlabIO.h"
#include "MatlabIO_WT.h"
#include "MatlabIO_getfield.h"
#include "myallocs.h"
#include "prnmsg.h"

#ifdef NoThreads

int pth_NRe_complex_1D_wavelet_rec(double *xrec, t_CWTvar *wdec, int L, int Tr, t_WaveletFamily *pWF) {

	NRe_complex_1D_wavelet_rec(xrec, wdec, L, Tr, pWF);
	return 0;
}

int pth_Ncomplex_1D_wavelet_rec(double complex *xrec, t_CWTvar *wdec, int L, int Tr, t_WaveletFamily *pWF) {

	Ncomplex_1D_wavelet_rec(xrec, wdec, L, Tr, pWF);
	return 0;
}

#else

#ifdef _MSC_VER
#include <windows.h>
#include <process.h>
#else
#include "pthread.h"
#endif

#define NTH 8

typedef	struct {
	double *xrec;
	t_CWTvar *wdec;
	t_WaveletFamily *pWF;
	int L, M;
} th_data_dbICWT_t;

void *thf_NRe_complex_1D_wavelet_rec(void *pvoid) {
	th_data_dbICWT_t *data = (th_data_dbICWT_t *)pvoid;

	NRe_complex_1D_wavelet_rec(data->xrec, data->wdec, data->L, data->M, data->pWF);

	return NULL;
}

int pth_NRe_complex_1D_wavelet_rec(double *xrec, t_CWTvar *wdec, int L, int Tr, t_WaveletFamily *pWF) {
	th_data_dbICWT_t data[NTH];
#ifdef _MSC_VER
	HANDLE thread_id[NTH];
#else
	pthread_t thread_id[NTH];
#endif
	int th, tr_ini, tr_end;

	tr_ini = 0;
	for (th=0; th<NTH; th++) {
		tr_end = (int)(ceil((th+1) * Tr/(double)NTH));
		if (tr_end > Tr) tr_end = Tr;
		data[th].xrec = xrec + L*tr_ini;
		data[th].wdec = wdec + tr_ini;
		data[th].pWF = pWF;
		data[th].L = L;
		data[th].M = tr_end-tr_ini;
		tr_ini = tr_end;
	}
#ifdef _MSC_VER
	for (th=0; th<NTH; th++) {
		thread_id[th] = (HANDLE)_beginthreadex(NULL, 0, (void *)thf_NRe_complex_1D_wavelet_rec, (void *)(&data[th]), 0, NULL);
		if (thread_id[th] == NULL) mexWarnMsgTxt("Create thread error");
	}
	for (th=0; th<NTH; th++) {
		WaitForSingleObject(thread_id[th], INFINITE);
		CloseHandle(thread_id[th]);
	}
#else
	for (th=0; th<NTH; th++) {
		if (pthread_create(&thread_id[th], NULL, (void *)thf_NRe_complex_1D_wavelet_rec, &data[th]) )
            mexWarnMsgTxt("Create thread error");
	}
	for (th=0; th<NTH; th++) pthread_join(thread_id[th], NULL);
#endif

	return 0;
}

void *thf_Ncomplex_1D_wavelet_rec(void *pvoid) {
	th_data_dbICWT_t *data = (th_data_dbICWT_t *)pvoid;

	Ncomplex_1D_wavelet_rec((double complex *)data->xrec, data->wdec, data->L, data->M, data->pWF);

	return NULL;
}

int pth_Ncomplex_1D_wavelet_rec(double complex *xrec, t_CWTvar *wdec, int L, int Tr, t_WaveletFamily *pWF) {
	th_data_dbICWT_t data[NTH];
#ifdef _MSC_VER
	HANDLE thread_id[NTH];
#else
	pthread_t thread_id[NTH];
#endif
	int th, tr_ini, tr_end;

	tr_ini = 0;
	for (th=0; th<NTH; th++) {
		tr_end = (int)(ceil((th+1) * Tr/(double)NTH));
		if (tr_end > Tr) tr_end = Tr;
		data[th].xrec = (double *)(xrec + L*tr_ini);
		data[th].wdec = wdec + tr_ini;
		data[th].pWF = pWF;
		data[th].L = L;
		data[th].M = tr_end-tr_ini;
		tr_ini = tr_end;
	}
#ifdef _MSC_VER
	for (th=0; th<NTH; th++) {
		thread_id[th] = (HANDLE)_beginthreadex(NULL, 0, (void *)thf_Ncomplex_1D_wavelet_rec, (void *)(&data[th]), 0, NULL);
		if (thread_id[th] == NULL) mexWarnMsgTxt("Create thread error");
	}
	for (th=0; th<NTH; th++) {
		WaitForSingleObject(thread_id[th], INFINITE);
		CloseHandle(thread_id[th]);
	}
#else
	for (th=0; th<NTH; th++) {
		if (pthread_create(&thread_id[th], NULL, (void *)thf_Ncomplex_1D_wavelet_rec, &data[th]) )
			mexWarnMsgTxt("Create thread error");
	}
	for (th=0; th<NTH; th++) pthread_join(thread_id[th], NULL);
#endif

	return 0;
}


#endif

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{
	t_WaveletDef WTdef;
	t_WaveletFamily *pWF;
	double  *pd1, *pd2;
	int ndim, cell_dims[2], m, L, M;
	mwSize dims[2];
	unsigned int s, i, ReOut;
	mxComplexity mxC;
	mxArray *pm;
		
	if (nrhs < 2) mexErrMsgTxt("Two inputs are required.");

	/* Read the wavelet variable, first part */

	if (prhs[0] == NULL) mexErrMsgTxt("wdec: No input data.");
	ndim = mxGetNumberOfDimensions(prhs[0]);

	if (ndim > 2) mexErrMsgTxt("wdec: Too many dimensions.");
	cell_dims[0] = mxGetM(prhs[0]);
	cell_dims[1] = mxGetN(prhs[0]);
	M = cell_dims[0]*cell_dims[1];

	/* Read the wavelet family inputs and build it. */

	WTdef.Family = GetIntegerField(prhs[1], "Family");
	mxC = (WTdef.Family > 0) ? mxREAL : mxCOMPLEX;
	WTdef.J = GetIntegerField(prhs[1], "J");
	WTdef.V = GetIntegerField(prhs[1], "V");
	WTdef.continuous = GetIntegerField(prhs[1], "continuous");
	ReOut = GetIntegerField(prhs[1], "RealInvOutput");
	if (WTdef.J < 1) mexErrMsgTxt("J: Must be a positive number.");
	if (WTdef.V < 1) mexErrMsgTxt("V: Must be a positive number.");
	if (WTdef.continuous != 0 && WTdef.continuous != 1) mexErrMsgTxt("continuous: Must be 0 or 1.");

	GetPositiveDoubleField (&WTdef.s0,  prhs[1], "s0");
	GetPositiveDoubleField (&WTdef.b0,  prhs[1], "b0");
	GetPositiveDoubleField (&WTdef.op1, prhs[1], "op1");

	if (nrhs >= 3) L = (int)mxGetScalar(prhs[2]);
	else {
		pd1 = mxGetPr(mxGetField(prhs[0], 0, "N"));
		L = (int)((int)*pd1 / WTdef.b0);
		if (WTdef.b0 >= 1) L = (int)*pd1;
	}
	WTdef.MWL = L;

	WTdef.convtype = GetIntegerFieldDefault (prhs[1], "convtype", 0);
	if (WTdef.convtype < 0 || WTdef.convtype > 1) mexErrMsgTxt("convtype: Method not defined.");
	
	/* Work a little bit. */
#ifndef NoThreads
	myallocs_init();
	prerror_init();
#endif

	pWF = CreateWaveletFamily (WTdef.Family, WTdef.J, WTdef.V, WTdef.MWL, WTdef.s0, WTdef.b0, WTdef.convtype, WTdef.op1, WTdef.continuous);
	if (pWF == NULL) mexErrMsgTxt("Wavelet Function not defined.");

	/* Build output variable. */
	if (M == 1) {dims[0] = 1; dims[1] = L;}
	else        {dims[0] = L; dims[1] = M;}
	if (ReOut == 1) {
		if (NULL == (plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL) ))
			mexErrMsgTxt("plhs[0]: Out of memory.");
	} else
		if (NULL == (plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxC) ))
			mexErrMsgTxt("plhs[0]: Out of memory.");

	if (mxC == mxREAL) {
		t_RWTvar *wdec;
		double *xrec;

		xrec = mxGetPr(plhs[0]);
		
		/* Read the wavelet variable, second part */
		wdec = CreateRealWaveletVarArray (pWF, L, M);
		for (m=0; m<M; m++) {
			pm = mxGetField(prhs[0], m, "d");
			for (s=0; s<wdec[m].S; s++) {
				pd1 = mxGetPr(mxGetCell(pm, s));
				pd2 = wdec[m].d[s];
				for (i=0; i<wdec[m].N[s]; i++) *pd2++ = *pd1++;
			}
		}

		Nreal_1D_wavelet_rec(xrec, wdec, L, M, pWF);
				
		DestroyRealWaveletVarArray(wdec, M);
	} else {
		t_CWTvar *wdec;
		double complex *xrec;

		/* Read the wavelet variable, second part */
		wdec = CreateComplexWaveletVarArray (pWF, L, M);
		for (m=0; m<M; m++) {
			pm = mxGetField(prhs[0], m, "d");
			for (s=0; s<wdec[m].S; s++)
				DA2CA_nomem(wdec[m].d[s], mxGetCell(pm, s));
		}
		if (ReOut == 0) {
			if (NULL == (xrec = (double complex *)mxMalloc(L*M*sizeof(double complex)) ))
				mexErrMsgTxt("xrec: Out of memory.");
			pth_Ncomplex_1D_wavelet_rec(xrec, wdec, L, M, pWF);
			CA2DA(plhs[0], xrec);
			mxFree(xrec);
		} else 
			pth_NRe_complex_1D_wavelet_rec(mxGetPr(plhs[0]), wdec, L, M, pWF);

		DestroyComplexWaveletVarArray(wdec, M);
	}
	
	plhs[1] = export_WaveletFamilyArray(&pWF, 1, 1);
	DestroyWaveletFamily(pWF);
#ifndef NoThreads
	prerror_destroy();
	myallocs_destroy();
#endif

}

#endif
