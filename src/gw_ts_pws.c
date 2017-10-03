#ifdef MATLAB
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "myallocs.h"
#include "prnmsg.h"
#include "MatlabIO_getfield.h"
#include "ts_pws1e_lib.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{
	const char *field_names[] = {"LS", "tsPWS", "LS_sim", "tsPWS_sim", "LS_misfit", "tsPWS_misfit"};
	t_tsPWS tsPWS = {-1, 0, 0, 4, 2., 1.0, PI*sqrt(2/log(2)), 2., 0., 0., 2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL};
	t_tsPWS_out out = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0};
	t_data      in;
	mxArray *pm;
	int ndim, ia1;
	unsigned int n, Tr, N;
	double *pd1, da1;
	float  *pf1;

	if (nrhs < 2) mexErrMsgTxt("Two inputs are required.");
	
	/*****************************/
	/* Read data related fields. */
	/*****************************/
	/* The Collection of traces. */
	pm = mxGetField(prhs[0], 0, "x");
	if (pm == NULL) mexErrMsgTxt("x: No input data.");
	ndim = mxGetNumberOfDimensions(pm);
	if (ndim > 2) mexErrMsgTxt("x: Too many dimensions.");
	if (!mxIsDouble(pm) || mxIsComplex(pm)) mexErrMsgTxt("x: Invalid type.");
	N  = mxGetM(pm);
	Tr = mxGetN(pm);
	pd1 = mxGetPr(pm);
	
	if (NULL == (in.sigall = (float *)mxMalloc(N*Tr*sizeof(float)) )) 
		mexErrMsgTxt("x: Out of memory.");
	pf1 = in.sigall;
	for (n=0; n<N*Tr; n++) pf1[n] = (float)pd1[n];
	
	/*  t_time field  */
	in.time = NULL;  /* Not used :-) */
	
	/* The trace header  */
	in.hdr.max = N;
	in.hdr.mtr = Tr;
	GetDoubleFieldDefault (&da1, prhs[0], "evla", 0); in.hdr.evla = da1;
	GetDoubleFieldDefault (&da1, prhs[0], "evlo", 0); in.hdr.evlo = da1;
	GetDoubleFieldDefault (&da1, prhs[0], "stla", 0); in.hdr.stla = da1;
	GetDoubleFieldDefault (&da1, prhs[0], "stlo", 0); in.hdr.stlo = da1;
	GetDoubleFieldDefault (&da1, prhs[0], "dt",   1); in.hdr.dt   = da1;
	GetDoubleFieldDefault (&da1, prhs[0], "beg",  0); in.hdr.beg  = da1;
	
	/* Reference trace */
	pm = mxGetField(prhs[0], 0, "reference");
	if (pm == NULL) in.reference = NULL;
	else {
		ndim = mxGetNumberOfDimensions(pm);
		if (ndim > 2) mexErrMsgTxt("reference: Too many dimensions.");
		if (!mxIsDouble(pm) || mxIsComplex(pm)) mexErrMsgTxt("reference: Invalid type.");
		n = mxGetM(pm)*mxGetN(pm);
		pd1 = mxGetPr(pm);
		if (n != N) {
			mexWarnMsgTxt("reference: Ignoring due to wrong length.");
			in.reference = NULL;
		} else {
			if (NULL == (in.reference = (float *)mxMalloc(N*sizeof(float)) )) 
				mexErrMsgTxt("reference: Out of memory.");
			pf1 = in.reference;
			for (n=0; n<N; n++) pf1[n] = (float)pd1[n];
		}
	}
	
	/*****************************/
	/* Configuration of tsPWS.   */
	/*****************************/
	/* Wavelet related fields    */
	pm = mxGetField(prhs[1], 0, "MexHat");
	if (pm != NULL) {
		tsPWS.type = -3;
		tsPWS.J    =  6;
		tsPWS.V    =  2;
		tsPWS.b0   =  0.5; 
		tsPWS.s0   =  1;
	}

	tsPWS.type = GetIntegerFieldDefault(prhs[1], "type", tsPWS.type);
	tsPWS.J = GetIntegerFieldDefault(prhs[1], "J", tsPWS.J);
	tsPWS.V = GetIntegerFieldDefault(prhs[1], "V", tsPWS.V);
	if (tsPWS.J < 1) mexErrMsgTxt("J: Must be a positive number.");
	if (tsPWS.V < 1) mexErrMsgTxt("V: Must be a positive number.");
	tsPWS.uni = GetIntegerFieldDefault(prhs[1], "uni", tsPWS.uni);
	if (tsPWS.uni != 0) tsPWS.uni = 1;
	GetPositiveDoubleFieldDefault (&tsPWS.s0,    prhs[1], "s0",    tsPWS.s0);
	GetPositiveDoubleFieldDefault (&tsPWS.b0,    prhs[1], "b0",    tsPWS.b0);
	GetPositiveDoubleFieldDefault (&tsPWS.Q,     prhs[1], "Q",     tsPWS.Q);
	
	GetPositiveDoubleFieldDefault (&tsPWS.cycle, prhs[1], "cycle", tsPWS.cycle);
	 
	GetPositiveDoubleFieldDefault (&tsPWS.w0,    prhs[1], "cyc",   tsPWS.w0);
	if ( mxGetField(prhs[1], 0, "cyc") ) tsPWS.w0 *= PI;
	GetPositiveDoubleFieldDefault (&tsPWS.w0,    prhs[1], "w0",    tsPWS.w0);
	
	if ( mxGetField(prhs[1], 0, "Q")      && tsPWS.w0set < 1) tsPWS.w0set = 1;
	if ( mxGetField(prhs[1], 0, "cycles") && tsPWS.w0set < 2) tsPWS.w0set = 2;
	if ( mxGetField(prhs[1], 0, "cyc")    && tsPWS.w0set < 3) tsPWS.w0set = 3;
	if ( mxGetField(prhs[1], 0, "w0")     && tsPWS.w0set < 4) tsPWS.w0set = 4;
	
	if ( mxGetField(prhs[1], 0, "V") )  tsPWS.lVfix  = 1;
	if ( mxGetField(prhs[1], 0, "b0") ) tsPWS.lb0fix = 1;
	if ( mxGetField(prhs[1], 0, "s0") ) tsPWS.ls0fix = 1;
	if ( mxGetField(prhs[1], 0, "verbose") ) tsPWS.verbose = 1;
	
	/* The other parameters */
	GetPositiveDoubleFieldDefault (&tsPWS.wu,    prhs[1], "wu", tsPWS.wu);
	tsPWS.lrm = GetIntegerFieldDefault(prhs[1], "rm", tsPWS.lrm);
	if (tsPWS.lrm) tsPWS.lrm = 1;
	tsPWS.fold = GetIntegerFieldDefault(prhs[1], "fold", tsPWS.fold);
	if (tsPWS.fold) tsPWS.fold = 1;
	tsPWS.unbiased = GetIntegerFieldDefault(prhs[1], "unbiased", tsPWS.unbiased);
	if (tsPWS.unbiased) tsPWS.unbiased = 1;
	tsPWS.convergence = GetIntegerFieldDefault(prhs[1], "convergence", tsPWS.convergence);
	if (tsPWS.convergence) tsPWS.convergence = 1;
	ia1 = GetIntegerFieldDefault(prhs[1], "Nmax", 0);
	tsPWS.Nmax = (ia1 < 0) ? 0 : ia1;
	ia1 = GetIntegerFieldDefault(prhs[1], "TwoStage", 0);
	tsPWS.Kmax = (ia1 < 0) ? 0 : ia1;

	/* Allocationg memory for results */
	plhs[0] = mxCreateStructMatrix(1, 1, sizeof(field_names)/sizeof(*field_names), field_names);
	mxSetField(plhs[0], 0, "LS",    mxCreateNumericMatrix(1, N, mxDOUBLE_CLASS, mxREAL));
	mxSetField(plhs[0], 0, "tsPWS", mxCreateNumericMatrix(1, N, mxDOUBLE_CLASS, mxREAL));
	if (tsPWS.convergence) {
		mxSetField(plhs[0], 0, "LS_sim",       mxCreateNumericMatrix(1, Tr, mxDOUBLE_CLASS, mxREAL));
		mxSetField(plhs[0], 0, "tsPWS_sim",    mxCreateNumericMatrix(1, Tr, mxDOUBLE_CLASS, mxREAL));
		mxSetField(plhs[0], 0, "LS_misfit",    mxCreateNumericMatrix(1, Tr, mxDOUBLE_CLASS, mxREAL));
		mxSetField(plhs[0], 0, "tsPWS_misfit", mxCreateNumericMatrix(1, Tr, mxDOUBLE_CLASS, mxREAL));
	} else {
		mxSetField(plhs[0], 0, "LS_sim",       mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL));
		mxSetField(plhs[0], 0, "tsPWS_sim",    mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL));
		mxSetField(plhs[0], 0, "LS_misfit",    mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL));
		mxSetField(plhs[0], 0, "tsPWS_misfit", mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL));
	}

	out.N     = N;
	out.mtr   = Tr;
	out.ls    = (float *)mxCalloc(N,sizeof(float));
	out.tsPWS = (float *)mxCalloc(N,sizeof(float));
	pm = mxGetField(plhs[0], 0, "LS_sim");       out.ls_sim       = mxGetPr(pm);
	pm = mxGetField(plhs[0], 0, "tsPWS_sim");    out.tsPWS_sim    = mxGetPr(pm);
	pm = mxGetField(plhs[0], 0, "LS_misfit");    out.ls_misfit    = mxGetPr(pm);
	pm = mxGetField(plhs[0], 0, "tsPWS_misfit"); out.tsPWS_misfit = mxGetPr(pm);


	/* Work a little bit. */
#ifndef NoThreads
	myallocs_init();
	prerror_init();
#endif
	tspws_main(&tsPWS, &out, &in);
#ifndef NoThreads
	prerror_destroy();
	myallocs_destroy();
#endif

	/* Cast float fields to double */
	pf1 = out.ls;
	pd1 = mxGetPr(mxGetField(plhs[0], 0, "LS"));
	for (n=0; n<N; n++) pd1[n] = pf1[n];
	pf1 = out.tsPWS;
	pd1 = mxGetPr(mxGetField(plhs[0], 0, "tsPWS"));
	for (n=0; n<N; n++) pd1[n] = pf1[n];
	
	mxFree(in.sigall);
	mxFree(in.reference);
	mxFree(out.ls);
	mxFree(out.tsPWS);
}

#endif
