#ifdef MATLAB
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "wavelet_v7.h"
#include "MatlabIO_WT.h"
#include "MatlabIO_getfield.h"
#include "myallocs.h"
#include "prnmsg.h"

#ifdef NoThreads

int pth_Ncomplex_1D_wavelet_dec(t_CWTvar *wdec, double *x, int L, int Tr, t_WaveletFamily *pWF) {

	Ncomplex_1D_wavelet_dec (wdec, x, L, Tr, pWF);
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
	t_CWTvar *wdec;
	double *x;
	t_WaveletFamily *pWF;
	int L, M;
} th_data_dbCWT_t;

void *thf_Ncomplex_1D_wavelet_dec(void *pvoid) {
	th_data_dbCWT_t *data = (th_data_dbCWT_t *)pvoid;

	Ncomplex_1D_wavelet_dec (data->wdec, data->x, data->L, data->M, data->pWF);

	return NULL;
}

int pth_Ncomplex_1D_wavelet_dec(t_CWTvar *wdec, double *x, int L, int Tr, t_WaveletFamily *pWF) {
	th_data_dbCWT_t data[NTH];
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
		data[th].wdec = wdec + tr_ini;
		data[th].x = x + L*tr_ini;
		data[th].pWF = pWF;
		data[th].L = L;
		data[th].M = tr_end-tr_ini;
		tr_ini = tr_end;
	}
#ifdef _MSC_VER
	for (th=0; th<NTH; th++) {
		thread_id[th] = (HANDLE)_beginthreadex(NULL, 0, (void *)thf_Ncomplex_1D_wavelet_dec, (void *)(&data[th]), 0, NULL);
		if (thread_id[th] == NULL) mexWarnMsgTxt("Create thread error");
	}
	for (th=0; th<NTH; th++) {
		WaitForSingleObject(thread_id[th], INFINITE);
		CloseHandle(thread_id[th]);
	}
#else
	for (th=0; th<NTH; th++) {
		if (pthread_create(&thread_id[th], NULL, (void *)thf_Ncomplex_1D_wavelet_dec, &data[th]))
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
	double *x;
	int ndim, dims[2], L, M;
	mxComplexity mxC;
		
	/* Read the input parameters. */
	if (nrhs < 2) mexErrMsgTxt("Two inputs are required.");
	
	if (prhs[0] == NULL) mexErrMsgTxt("x: No input data.");
	ndim = mxGetNumberOfDimensions(prhs[0]);
	if (ndim > 2) mexErrMsgTxt("x: Too many dimensions.");
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) mexErrMsgTxt("x: Invalid type.");
	dims[0] = mxGetM(prhs[0]);
	dims[1] = mxGetN(prhs[0]);
	if (dims[0] == 1) {L = dims[1]; M = 1;}
	else              {L = dims[0]; M = dims[1];}
	x = mxGetPr(prhs[0]);

	WTdef.Family = GetIntegerField(prhs[1], "Family");
	mxC = (WTdef.Family > 0) ? mxREAL : mxCOMPLEX;
	WTdef.J = GetIntegerField(prhs[1], "J");
	WTdef.V = GetIntegerField(prhs[1], "V");
	WTdef.continuous = GetIntegerField(prhs[1], "continuous");
	WTdef.MWL = L;
	if (WTdef.J < 1) mexErrMsgTxt("J: Must be a positive number.");
	if (WTdef.V < 1) mexErrMsgTxt("V: Must be a positive number.");
	if (WTdef.continuous != 0 && WTdef.continuous != 1) mexErrMsgTxt("continuous: Must be 0 or 1.");

	GetPositiveDoubleField (&WTdef.s0,  prhs[1], "s0");
	GetPositiveDoubleField (&WTdef.b0,  prhs[1], "b0");
	GetPositiveDoubleField (&WTdef.op1, prhs[1], "op1");

	WTdef.convtype = GetIntegerFieldDefault (prhs[1], "convtype", 0);
	if (WTdef.convtype < 0 || WTdef.convtype > 1) mexErrMsgTxt("convtype: Method not defined.");
	
	/* Work a little bit. */
#ifndef NoThreads
	myallocs_init();
	prerror_init();
#endif
	pWF = CreateWaveletFamily (WTdef.Family, WTdef.J, WTdef.V, WTdef.MWL, WTdef.s0, WTdef.b0, WTdef.convtype, WTdef.op1, WTdef.continuous);
	if (pWF == NULL) mexErrMsgTxt("Wavelet Function not defined.");

	if (mxC == mxREAL) {
		t_RWTvar *wdec;

	    wdec = CreateRealWaveletVarArray (pWF, L, M);

		Nreal_1D_wavelet_dec (wdec, x, L, M, pWF);
		plhs[0] = export_RWT(wdec, 1, M);
		
		DestroyRealWaveletVarArray(wdec, M);
	} else {
		t_CWTvar *wdec;
		wdec = CreateComplexWaveletVarArray (pWF, L, M);

		/* The work */
		pth_Ncomplex_1D_wavelet_dec (wdec, x, L, M, pWF);
		plhs[0] = export_CWT(wdec, 1, M);

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
