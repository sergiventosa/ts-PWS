#ifndef SAC2BIN_H
#define SAC2BIN_H

#include <stdint.h>

typedef struct {
	char     method[8]; /* Cross-correlation method. */
	char     net1[8];   /* Net 1 name.               */
	char     sta1[8];   /* Station 1 name.           */
	char     loc1[8];   /* Location/hole 1 id.       */
	char     chn1[8];   /* Channel id.               */
	char     net2[8];   /* Net 1 name.               */
	char     sta2[8];   /* Station 2 name.           */
	char     loc2[8];   /* Location/hole 1 id.       */
	char     chn2[8];   /* Channel id.               */
	float    stlat1;    /* Latitud of st. 1.         */
	float    stlon1;    /* Longitud of st. 1.        */
	float    stel1;     /* Elevetion of st. 1.       */
	float    stlat2;    /* Latitud of st. 2.         */
	float    stlon2;    /* Longitud of st. 2.        */
	float    stel2;     /* Elevetion of st. 2.       */
	uint32_t nlags;     /* Number of lags.           */
	uint32_t nseq;      /* Number of sequences.      */
	float    tlength;   /* Length of the corr. seq. (s) */
	float    lag1;      /* Lowest lag time (s).      */
	float    lag2;      /* Highest lag time (s).     */
} t_ccheader;

#endif
