#ifdef MATLAB
#include "MatlabIO.h"
#include "MatlabIO_getfield.h"

int GetIntegerField (const mxArray *p, const char *fieldname) {
	mxArray *pm;

	pm = mxGetField(p, 0, fieldname);
	if (pm == NULL) mexerror(fieldname,"Nonexistent field.");
	if (!mxIsDouble(pm) || mxIsComplex(pm) || mxGetN(pm)*mxGetM(pm) != 1) 
		mexerror(fieldname,"Input must be a scalar.");
	return (int)mxGetScalar(pm);
}

unsigned int GetUIntegerField (const mxArray *p, const char *fieldname) {
	mxArray *pm;

	pm = mxGetField(p, 0, fieldname);
	if (pm == NULL) mexerror(fieldname,"Nonexistent field.");
	if (!mxIsDouble(pm) || mxIsComplex(pm) || mxGetN(pm)*mxGetM(pm) != 1) 
		mexerror(fieldname,"Input must be a scalar.");
	return (unsigned int)mxGetScalar(pm);
}

int GetIntegerFieldDefault (const mxArray *p, const char *fieldname, int def) {
	mxArray *pm;

	pm = mxGetField(p, 0, fieldname);
	if (pm == NULL) return def;
	else {
		if (!mxIsDouble(pm) || mxIsComplex(pm) || mxGetN(pm)*mxGetM(pm) != 1) 
			mexerror(fieldname,"Input must be a scalar.");
		return (int)mxGetScalar(pm);
	}
}

void GetDoubleField (double *out, const mxArray *p, const char *fieldname) {
	mxArray *pm;

	pm = mxGetField(p, 0, fieldname);
	if (pm == NULL) mexerror(fieldname,"Nonexistent field.");
	if (!mxIsDouble(pm) || mxIsComplex(pm) || mxGetN(pm)*mxGetM(pm) != 1) 
		mexerror(fieldname,"Input must be a scalar.");
	*out = mxGetScalar(pm);
}

void GetPositiveDoubleField (double *out, const mxArray *p, const char *fieldname) {
	mxArray *pm;

	pm = mxGetField(p, 0, fieldname);
	if (pm == NULL) mexerror(fieldname,"Nonexistent field.");
	if (!mxIsDouble(pm) || mxIsComplex(pm) || mxGetN(pm)*mxGetM(pm) != 1) 
		mexerror(fieldname,"Input must be a scalar.");
	*out = mxGetScalar(pm);
	if (*out < 0) mexerror(fieldname,"Must be a positive number.");
}

void GetPositiveDoubleFieldDefault (double *out, const mxArray *p, const char *fieldname, double def) {
	mxArray *pm;

	pm = mxGetField(p, 0, fieldname);
	if (pm == NULL) *out = def;
	else {
		if (!mxIsDouble(pm) || mxIsComplex(pm) || mxGetN(pm)*mxGetM(pm) != 1) 
			mexerror(fieldname,"Input must be a scalar.");
		*out = mxGetScalar(pm);
		if (*out < 0) mexerror(fieldname,"Must be a positive number.");
	}
}

void GetDoubleFieldDefault (double *out, const mxArray *p, const char *fieldname, double def) {
	mxArray *pm;

	pm = mxGetField(p, 0, fieldname);
	if (pm == NULL) *out = def;
	else {
		if (!mxIsDouble(pm) || mxIsComplex(pm) || mxGetN(pm)*mxGetM(pm) != 1) 
			mexerror(fieldname,"Input must be a scalar.");
		*out = mxGetScalar(pm);
	}
}

#endif
