#ifndef MATLABIOGF_H
#define MATLABIOGF_H

#ifdef MATLAB
#include "mex.h"
#include "matrix.h"

int GetIntegerField (const mxArray *p, const char *fieldname);
unsigned int GetUIntegerField (const mxArray *p, const char *fieldname);
int GetIntegerFieldDefault (const mxArray *p, const char *fieldname, int def);
void GetDoubleField (double *out, const mxArray *p, const char *fieldname);
void GetPositiveDoubleField (double *out, const mxArray *p, const char *fieldname);
void GetPositiveDoubleFieldDefault (double *out, const mxArray *p, const char *fieldname, double def);
void GetDoubleFieldDefault (double *out, const mxArray *p, const char *fieldname, double def);

#endif

#endif

