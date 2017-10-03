#ifdef MATLAB
#include <string.h>

#include "matrix.h"
#include "MatlabIO.h"

void mexerror(const char *str1, const char *str2) {
	char *s;
	unsigned int l;

	l = strlen(str1) + strlen(str2) + 3;
	s = (char *)mxMalloc(l*sizeof(char));
	if (s == NULL) mexErrMsgTxt("mexerror: Out of memory");
	sprintf(s, "%s: %s", str1, str2);
	mexErrMsgTxt(s);
}

complex *DA2CA(const mxArray *mx) {
	double complex *pc, *pcaux;
	double *xr, *xi;
	int n, N;

	N = mxGetNumberOfElements(mx);
	if (mxIsComplex(mx)) {
		xr = mxGetPr(mx);
		xi = mxGetPi(mx);
	} else return(NULL);
	pc = (complex *)mxMalloc(N*sizeof(complex));
	if (pc != NULL) {
		pcaux = pc;
		for (n=0; n<N; n++)
			pcaux[n] = xr[n] + I * xi[n];
	}
	return(pc);
}

void DA2CA_nomem(complex *pc, const mxArray *mx) {
	double *xr, *xi;
	int n, N;

	N = mxGetNumberOfElements(mx);
	if (mxIsComplex(mx)) {
		xr = mxGetPr(mx);
		xi = mxGetPi(mx);
		if (pc != NULL) {
			for (n=0; n<N; n++)
				pc[n] = xr[n] + I * xi[n];
		}
	} else {
		xr = mxGetPr(mx);
		if (pc != NULL) {
			for (n=0; n<N; n++)
				pc[n] = xr[n];
		}
	}
}

int CA2DA(mxArray *mx, complex *xc) {
	double *xr, *xi;
	int n, N;

	if (mx != NULL) {
		xr = mxGetPr(mx);
		xi = mxGetPi(mx);
	} else return 1;
	N = mxGetNumberOfElements(mx);
	for (n=0; n<N; n++) {
		xr[n] = creal(xc[n]);
		xi[n] = cimag(xc[n]);
	}
	return 0;
}

void check_double_scalar(const mxArray *p, char *msg) {
	if (!mxIsDouble(p) || mxIsComplex(p) || mxGetN(p)*mxGetM(p) != 1) 
		mexErrMsgTxt(msg);
}

#endif
