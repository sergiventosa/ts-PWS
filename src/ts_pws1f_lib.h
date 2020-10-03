/*****************************************************************************/
/* Program to perform time-scale phase-weiged stack (ts-pws) using frames    */
/* of continuous wavelets.                                                   */
/*                                                                           */
/* Author: Sergi Ventosa Rahuet (sergiventosa@hotmail.com)                   */
/*                                                                           */
/* The actual library computing ts_pws & friends.                            */
/*****************************************************************************/
#ifndef TSPWS_LIB
#define TSPWS_LIB

#include <stdio.h>
#include <time.h>

#ifndef PI
#define	PI 3.14159265358979328
#endif

typedef struct {
	int           type;       /* Wavelet Family Id.                           */
	unsigned int  uni;        /* Dyadic (0, default) / Uniform (1) sampling   */
	unsigned int  J;          /* Number of scales (default 6).                */
	unsigned int  V;          /* Number of voices per octave (default 4).     */
	double        s0;         /* First scale (default 2)                      */
	double        b0;         /* Sampling Period at scale = 1 (default 1)     */
	double        w0;         /* Mother-wavelet central frequency             */
	                          /* (default pi*sqrt(2/log(2)).                  */
	double        wu;         /* (TS-)PWS power (default 2)                   */
	double        fmin;       /* Minimum frequency analyzed.                  */
	double        Q;          /* Quality.                                     */
	double        cycle;      /* Number of cycles.                            */
	int           w0set;      /* (Internal) 0 not set, 1 Q, 2 cycle & 3 w0.   */
	int           lrm;        /* Remove the mean (default false).             */
	int           bin;        /* Input in binary format (instead of sacs).    */
	int           lkinst;
	int           lVfix;
	int           ls0fix;
	int           lb0fix;
	int           verbose;
	int           fold;
	int           unbiased;
	int           convergence;
	unsigned int  subsmpl_N;
	double        subsmpl_p;
	unsigned int  jackknife_n;
	unsigned int  jackknife_d;
	unsigned int  obin;
	int           AllSteps;   /* When AllSteps & convergence are true, one    */
	                          /* output is written for each stacking steps.   */
	unsigned int  Nmax;       /* Using Nmax sequence only.                    */
	unsigned int  Kmax;       /* Two-stage with Kmax sequence in the ts-PWS.  */
	char          *kinst;     /* String to be written on the instrument field.*/
	char          *filein;    /* Name of the file containing a list of sac    */
	                          /* filenames or the sequences in binary format. */
	char          *fileout;   /* Common part of the output file names.        */
	char          *fileconv;  /* Alternative reference trace for convergence. */
} t_tsPWS;

typedef struct {
	int          max;     /* Number of samples.       */
	unsigned int mtr;     /* Number of sequences.     */
	float        evla;    /* Source (sta1) latitud.   */
	float        evlo;    /* Source longitud.         */
	float        stla;    /* Receiver (sta2) latitud. */
	float        stlo;    /* Receiver longitud.       */
	float        stel;    /* Receiver elevation.      */
	float        dt;      /* Sampling period.         */
	float        beg;     /* Starting time.           */
	char         net1[9]; /* Net name of sta1         */
	char         sta1[9]; /* Station 1 name.          */
	char         loc1[9]; /* Location/hole 1 id.      */
	char         chn1[9]; /* Channel id of sta1.      */
	char         net2[9]; /* Net name of sta2.        */
	char         sta2[9]; /* Station 2 name.          */
	char         loc2[9]; /* Location/hole 1 id.      */
	char         chn2[9]; /* Channel id of sta2.      */
} t_hdr;

typedef struct {
	float        *ls;             /* Linear stack.                                 */
	float        *tsPWS;          /* ts-PWS.                                       */
	double       *ls_sim;         /* Convergence (similarity) of the linear stack. */
	double       *tsPWS_sim;      /* Convergence (similarity) of the ts-PWS.       */
	double       *ls_misfit;      /* Convergence (misfit) of the linear stack.     */
	double       *tsPWS_misfit;   /* Convergence (misfit) of the ts-PWS.           */
	float        *ls_steps;
	float        *tsPWS_steps;
	float        **ls_subsmpl;    /* Subsampled linear stacks.                     */ 
	float        **tsPWS_subsmpl; /* Subsampled ts-PWS.                            */ 
	unsigned int *mtr_subsmpl;    /* Number of sequences involved in the subsmple. */
	unsigned int M;               /* Number of subsample realizations (subsmpl_M)  */
	unsigned int N;               /* Length of the above sequences.                */
	unsigned int mtr;             /* Number of sequences.                          */
} t_tsPWS_out;

typedef struct {
	float  *sigall;    /* Time sequences.                              */
	time_t *time;      /* Starting times.                              */
	float  *lag0;      /* Lag time due to trace misaligment.           */
	t_hdr   hdr;       /* Header.                                      */
	float  *reference; /* Alternative reference trace for convergence. */
} t_data;

int tspws_main(t_tsPWS *tspws, t_tsPWS_out *out, t_data *in);

#endif
