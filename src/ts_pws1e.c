/*****************************************************************************/
/* Program to perform time-scale phase-weiged stack (ts-PWS) using frames    */
/* of continuous wavelets, including the two-stage stack and the unbiased    */
/* phase coherence strategies. (Ventosa et al. (GJI, 2017)                   */
/*                                                                           */
/*  Ventosa, S., Schimmel, M., & Stutzmann, E., 2017. Extracting surface     */
/*  waves, hum and normal modes: Time-scale phase-weighted stack and beyond, */
/*  Geophysical Journal International, doi: 10.1093/gji/ggx284               */
/*                                                                           */
/* Author: Sergi Ventosa Rahuet (sergiventosa(at)hotmail.com)                */
/*                                                                           */
/* Main features:                                                            */
/* - Do the ts-pws using the Morlet wavelet or the Mexican hat wavelets      */
/*   implement using frames of wavelets for the time-frequency expansion.    */
/* - All parameters are automatically set, when w0, Q or number of cycles is */
/*   used the parameters descriving the frame are changed accordindly.       */
/* - OpenMP is used in the forward wavelet transform.                        */
/* **** 2012 ****                                                            */
/* Jul13-17: First version based on time-frequency phase-weighed stack       */
/*           previously developed by Martin Schimmel.                        */
/* **** 2016 ****                                                            */
/* Jan20 (1a) Migration to C99 double complex types.                         */
/* May09 (1b) Change to dynamic memory.                                      */
/* Jun08 (1c) Minor corrections and additions to the parameters.             */
/* Nov03 (1d) New features, and minor corrections.                           */
/*            - Support for a binary file format containing all the          */
/*              sequences to be stacked.                                     */
/*            - Support for Complex Mexican Hat wavelet.                     */
/*            - New parameters, such as folding and file output names.       */
/* **** 2017 ****                                                            */
/* Feb01 (1e) New features, and minor corrections.                           */
/*            - Convergence figures for the linear stack and ts-PWS.         */
/*            - Two-stage stacking.                                          */
/*            - Unbiased phase coherence.                                    */
/* Sep20 (1e) Documentation.                                                 */
/* **** 2018 ****                                                            */
/* Jan17 (1e) Verbose now shows the output file names for the linear stack   */
/*            and the ts-PWS.                                                */ 
/* Feb15 (1e) Provide more information errors & warnings.                    */
/*            Undefined event and station locations does not stop execution. */
/*            The fold parameter is ignored on non-symmetric sequences.      */
/* Mar09 (1e) Bug corrections:                                               */
/*            - Corrected the errors and warnings produced when using the    */
/*            default compiler in MAC.                                       */
/*            - The fold parameter is now ignored on non-symmetric sequences.*/
/*            - Odd-length sequences are now well folded (the zero-lag       */
/*            sample was halved).                                            */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sacio.h>
#include <math.h>

#include "sac2bin.h"
#include "myallocs.h"
#include "ts_pws1e_lib.h"

void infooo();
void usage();

int ReadData (float **sigall, time_t **time, t_hdr *hdr, char *filein, int bin, int verbose);
int wrsac(char *filename, float *y, t_hdr *hdr, char *kstnm, float user0);

void DestroyFilelist(char *p[]);
int CreateFilelist (char **filename[], unsigned int *Tr, char *filelist);

int isuint (const char *str);
int isdouble (const char *str);

int RDuint   (unsigned int* const x, const char *str);
int RDdouble (double* const x, const char *str);

int error_header (char *filename, int nerr) {
	printf ("\a tspws_main: Error reading %s header (nerr=%d).\n", filename, nerr); 
	return 2;
}

int error_missing_info (char *filename, char *nick, int nerr) {
	printf ("\a tspws_main: Error reading %s header (nerr=%d), %s is not defined!\n", filename, nerr, nick); 
	return 2;
}

void warning_missing_info (char *filename, char *nick, int verbose) {
	if (verbose)
		printf ("\a tspws_main: %s is not defined in %s.\n", filename, nick); 
}

void strcat3 (char *sout, const char *s1, const char *s2, const char *s3);

/********************************************************************************************/
/* Main function: Reads the parameters and the data, allocates memory and save the results. */
/********************************************************************************************/
int main(int argc, char *argv[]) {
	t_tsPWS tspws = {-1, 0, 0, 4, 2., 1.0, PI*sqrt(2/log(2)), 2., 0., 0., 2., 
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL};  /* Default parameters. */
	t_tsPWS_out out;
	t_data in;
	int i, er;
	char *filename=NULL;

	/********************/
	/* Read parameters. */
	/********************/
	if (argc == 1) {
		usage();
		return 0;
	}
	if (argc > 1) {
		tspws.filein = argv[1];
		if (!strncmp(tspws.filein, "info", 4)) {
			infooo();
			return 0;
		}
	}
	
	for (i=2; i<argc; i++) {
		if (!strncmp(argv[i], "wu=", 3))        er = RDdouble(&tspws.wu,   argv[i] + 3);
		else if (!strncmp(argv[i], "J=",    2)) er = RDuint(&tspws.J,      argv[i] + 2);
		else if (!strncmp(argv[i], "Nmax=", 5)) er = RDuint(&tspws.Nmax,   argv[i] + 5);
		else if (!strncmp(argv[i], "fmin=", 5)) er = RDdouble(&tspws.fmin, argv[i] + 5);
		else if (!strncmp(argv[i], "V=",  2)) { er = RDuint(&tspws.V,      argv[i] + 2); tspws.lVfix  = 1; }
		else if (!strncmp(argv[i], "s0=", 3)) { er = RDdouble(&tspws.s0,   argv[i] + 3); tspws.ls0fix = 1; }
		else if (!strncmp(argv[i], "b0=", 3)) { er = RDdouble(&tspws.b0,   argv[i] + 3); tspws.lb0fix = 1; }
		else if (!strncmp(argv[i], "uni", 3))      tspws.uni      =  1; /* True */
		else if (!strncmp(argv[i], "rm",  2))      tspws.lrm      =  1; /* True */
		else if (!strncmp(argv[i], "bin",          3)) tspws.bin         =  1; /* True */
		else if (!strncmp(argv[i], "verbose",      7)) tspws.verbose     =  1; /* True */
		else if (!strncmp(argv[i], "fold",         4)) tspws.fold        =  1; /* True */
		else if (!strncmp(argv[i], "unbiased",     8)) tspws.unbiased    =  1; /* True */
		else if (!strncmp(argv[i], "MexHat",       6)) tspws.type        = -3;
		else if (!strncmp(argv[i], "osac=",        5)) tspws.fileout     = argv[i] + 5;
		else if (!strncmp(argv[i], "AllSteps",     8)) tspws.AllSteps    =  1; /* True */
		else if (!strncmp(argv[i], "TwoStage", 8)) {
			tspws.Kmax =  10; /* The default value. */
			if (!strncmp(argv[i], "TwoStage=", 9)) er = RDuint(&tspws.Kmax, argv[i] + 9);
		} else if (!strncmp(argv[i], "convergence", 11)) {
			tspws.convergence =  1; /* True */
			if (!strncmp(argv[i], "convergence=", 12)) tspws.fileconv    = argv[i] + 12;
		} else if (!strncmp(argv[i], "Q=",  2)) {
			er = RDdouble(&tspws.Q, argv[i] + 2);
			if (tspws.w0set < 1) tspws.w0set = 1;
		} else if (!strncmp(argv[i], "cycles=", 7)) {
			er = RDdouble(&tspws.cycle, argv[i] + 7);
			if (tspws.w0set < 2) tspws.w0set = 2;
		} else if (!strncmp(argv[i], "cyc=", 4)) {  /* Because legacy with Martin's code. */
			if (tspws.w0set < 3) {
				tspws.w0set = 3;
				er = RDdouble(&tspws.w0, argv[i] + 4);
				tspws.w0 *= PI;
			}
		} else if (!strncmp(argv[i], "w0=", 3)) {
			er = RDdouble(&tspws.w0, argv[i] + 3);
			if (tspws.w0set < 4) tspws.w0set = 4;
		} else if (!strncmp(argv[i], "kstnm=", 6)) {
			tspws.kstnm  = argv[i] + 6;
			tspws.lkstnm = 1; /* True */
		} else if (!strncmp(argv[i], "info", 4)) {
			infooo();
			return 0;
		}
	}
	
	/********************/
	/* Initializations. */
	/********************/
	/* Read data and header */
	memset(&in, 0, sizeof(in));
	er = ReadData (&in.sigall, &in.time, &in.hdr, tspws.filein, tspws.bin, tspws.verbose);
	if (tspws.convergence && tspws.fileconv) {
		if (NULL == (in.reference = (float *)mycalloc(in.hdr.max, sizeof(float)) )) er = 4;
	} else in.reference = NULL;
	if (er) return er;
	
	/* Allocate memory */
	memset(&out, 0, sizeof(out));
	out.N   = in.hdr.max;
	out.mtr = (tspws.Nmax) ? tspws.Nmax : in.hdr.mtr;
	if (NULL == (out.ls    = (float *)mycalloc(in.hdr.max, sizeof(float)) )) er = 4;
	if (NULL == (out.tsPWS = (float *)mycalloc(in.hdr.max, sizeof(float)) )) er = 4;
	if (tspws.convergence) {
		if (NULL == (out.ls_sim       = (double *)mymalloc(in.hdr.mtr * sizeof(double)) )) er = 4;
		if (NULL == (out.tsPWS_sim    = (double *)mymalloc(in.hdr.mtr * sizeof(double)) )) er = 4;
		if (NULL == (out.ls_misfit    = (double *)mymalloc(in.hdr.mtr * sizeof(double)) )) er = 4;
		if (NULL == (out.tsPWS_misfit = (double *)mymalloc(in.hdr.mtr * sizeof(double)) )) er = 4;
		if (tspws.AllSteps) {
			if (NULL == (out.ls_steps    = (float *)mymalloc(in.hdr.mtr*in.hdr.max * sizeof(float)) )) er = 4;
			if (NULL == (out.tsPWS_steps = (float *)mymalloc(in.hdr.mtr*in.hdr.max * sizeof(float)) )) er = 4;
		}
	}
	if (er == 4) {
		myfree(out.ls);
		myfree(out.tsPWS);
		myfree(out.ls_sim);
		myfree(out.tsPWS_sim);
		myfree(out.ls_misfit);
		myfree(out.tsPWS_misfit);
		myfree(out.ls_steps);
		myfree(out.tsPWS_steps);
		printf ("main: Out of memory when reading %s (npts = %d)\n", tspws.filein, in.hdr.max);
		return er; 
	}
	
	if (tspws.fileconv) {
		t_hdr hdr2;
		int ia1;
		
		filename = tspws.fileconv;
		rsach (filename, &er, strlen(filename));
		if (!er) {
			sac_warning_off ();
			getnhv ("npts",  &hdr2.max,  &er, strlen("npts"));
		}
		if (!er) {
			rsac1(filename, in.reference, &ia1, &hdr2.beg, &hdr2.dt, &hdr2.max, &er, strlen(filename));
			if (in.hdr.max != ia1) {
				printf("The reference for convergence has a different length (%d:%d)\n", in.hdr.max, ia1);
				er = 1;
			}
		}
		if (er) { /* In case of any error, ignore the reference trace. */
			myfree(in.reference);
			in.reference   = NULL;
			tspws.fileconv = NULL;
			er = 0;
		}
	}
	
	/*****************/
	/* The main job. */
	/*****************/
	er = tspws_main(&tspws, &out, &in); /* The one who make the job. */
	
	/***********************/
	/* Save, clean & quit. */
	/***********************/
	/* Save */
	if (!er) {
		t_hdr hdr;
		unsigned int n, ua1;
		float *sig;
		FILE *fid;
		
		memcpy(&hdr, &in.hdr, sizeof(t_hdr));
		
		/* Write sac files */
		n = 0;
		if (tspws.fold) {
			if (2*hdr.beg + (hdr.max-1)*hdr.dt < 0.5*hdr.dt) {
				n = hdr.max/2;            /* First sample to be saved. */
				hdr.max = (hdr.max+1)/2;  /* Number of samples   "   . */
				hdr.beg += hdr.dt*n;      /* Update beg time.          */
			} else printf("Warning: Folding ignored. B = %f, E = %f\n", hdr.beg, hdr.beg + (hdr.max-1)*hdr.dt);
		}
		
		ua1 = 25;
		if (tspws.fileout) ua1 += strlen(tspws.fileout);
		filename = (char *)mymalloc(ua1*sizeof(char));
		
		/* Linear stack */
		if (tspws.fileout) strcat3 (filename, "tl_", tspws.fileout, ".sac");
		else strcpy(filename, "tl.sac");
		
		if (tspws.verbose) {
			printf("  Output files:\n");
			printf("     Linear stack: %s\n", filename);
		}
		
		sig = out.ls + n;
		if (tspws.lkstnm)
			wrsac(filename, sig, &hdr, tspws.kstnm, (float)out.mtr);
		else wrsac(filename, sig, &hdr, "t-lin", (float)out.mtr);
		
		/* TS-PWS */
		if (tspws.fileout) strcat3 (filename, "ts_pws_", tspws.fileout, ".sac");
		else strcpy(filename, "ts_pws.sac");
		
		sig = out.tsPWS + n;
		if (tspws.lkstnm)
			wrsac(filename, sig, &hdr, tspws.kstnm, (float)out.mtr);
		else wrsac(filename, sig, &hdr, "ts_pws", (float)out.mtr);
		
		if (tspws.verbose) 
			printf("     ts-PWS:       %s\n", filename);
		
		/* Convergence binary files. */
		if (tspws.convergence) {
			/* Convergence on similarity for the tsPWS. */
			strcat3 (filename, "ts_pws_", tspws.fileout, "_convergence");
			fid = fopen(filename, "w");
			fwrite (out.tsPWS_sim, sizeof(double), out.mtr, fid);
			fclose (fid);
			/* Convergence on similarity for the linear stack. */
			strcat3 (filename, "tl_", tspws.fileout, "_convergence");
			fid = fopen(filename, "w");
			fwrite (out.ls_sim, sizeof(double), out.mtr, fid);
			fclose (fid);
			/* Convergence on misfit for the tsPWS. */
			strcat3 (filename, "ts_pws_", tspws.fileout, "_misfit");
			fid = fopen(filename, "w");
			fwrite (out.tsPWS_misfit, sizeof(double), out.mtr, fid);
			fclose (fid);
			/* Convergence on misfit for the linear stack. */
			strcat3 (filename, "tl_", tspws.fileout, "_misfit");
			fid = fopen(filename, "w");
			fwrite (out.ls_misfit, sizeof(double), out.mtr, fid);
			fclose (fid);
			/* Write all the stacking steps. */
			if (out.tsPWS_steps) {
				strcat3 (filename, "ts_pws_", tspws.fileout, "_steps");
				fid = fopen(filename, "w");
				fwrite (out.tsPWS_steps, sizeof(float), out.mtr*out.N, fid);
				fclose (fid);
			}
			if (out.ls_steps) {
				strcat3 (filename, "tl_", tspws.fileout, "_steps");
				fid = fopen(filename, "w");
				fwrite (out.ls_steps, sizeof(float), out.mtr*out.N, fid);
				fclose (fid);
			}
		}
		myfree(filename);
	}
	/* clean */
	myfree(out.ls);
	myfree(out.tsPWS);
	myfree(out.ls_sim);
	myfree(out.tsPWS_sim);
	myfree(out.ls_misfit);
	myfree(out.tsPWS_misfit);
	myfree(out.ls_steps);
	myfree(out.tsPWS_steps);
	myfree(in.sigall);
	myfree(in.time);
	myfree(in.reference);
	
	return 0;
}

void infooo() {
	puts("This program performs a PWS (Schimmel & Paulssen, 1997) in the time-scale domain employing a Continuous Wavelet Transfrom (Ventosa et al., 2017). This is similar in spirit to tf_pws.f which employs the S-Transform (Schimmel & Gallart, 2007) but with a much lower CPU and memory requirements. The additional option of TwoStage and unbiased may improve results.");
	puts("Information on the relation between Wavelet and S transfrons is published in Ventosa et al. (2008).\n");
	puts("Schimmel M., & H. Paulssen, 1997. Noise reduction and detection of weak, coherent signals through phase weighted stacks, Geophysical Journal International, 130, 497-505. doi:10.1111/j.1365-246X.1997.tb05664.x");
	puts("Schimmel M., & J. Gallart, 2007. Frequency-dependent phase coherence for noise suppression in seismic array data, Journal of Geophysical Research, 112, B04303, doi:10.1029/2006JB004680");
	puts("Ventosa et al., 2008. S-transform from a wavelets point of view, IEEE Transactions on Signal Processing, 56, 2771-2780, doi:10.1109/TSP.2008.917029");
	puts("Ventosa, S., Schimmel, M., & Stutzmann, E., 2017. Extracting surface waves, hum and normal modes: Time-scale phase-weighted stack and beyond, Geophysical Journal International, 211, 30-44, doi:10.1093/gji/ggx284"); 
	puts("AUTHOR: Sergi Ventosa Rahuet (sergiventosa(at)hotmail.com)");
	puts("Last modification: 09/03/2018\n");
}

void usage() {
	puts("\nUSAGE: ts_pws filename parameter_list");
	puts("  filename: One file name (SAC) per line, traces must have same begin (b),");
	puts("            number of samples (nsmpl) and sampling interval (dt).");
	puts("  parameters are optional and can be provided in arbitrary order without any");
	puts("  blanck around = .");
	puts("");
	puts("Most commonly used parameters");
	puts("  cyc=   : Number of cycles defined at the one sigma (from -sigma to +sigma).");
	puts("           This is equivalent to w0=cyc*pi in the (default) Morlet wavelet.");
	puts("  fmin=  : Minimum frequency analyzed (Hz).");
	puts("  info   : write background and main references to screen.");
	puts("           Just type: ts_pws info.");
	puts("  kstnm= : optionally specify sac-header (char*8).");
	puts("  rm     : remove the mean from input data.");
	puts("  verbose: Print the parameters used and files stacked.");
	puts("  wu=    : PWS power.");
	puts("");
	puts("Parameters controlling the discretization of the Continuous Wavelet Transform.");
	puts("  b0=    : Sampling at the scale 1.");
	puts("  J=     : Number of octaves. Overwrites fmin.");
	puts("  s0=    : Lowest scale sampled. Overwrites fmin.");
	puts("  uni    : Use uniform sampling instead of the default dyadic sampling.");
	puts("  V=     : Number of voices per octaves.");
	puts("");
	puts("Parameters controlling the mother wavelet/function shape.");
	puts("  Chose one of w0, cyc, cycles or Q to control the time-frequency resolution.");
	puts("  When w0 (equivalently Q, cycles or cyc) is changed, V and b0 are adapted accordingly.");
	puts("  w0=    : Central frequency of the Morlet wavelet. Overwrites cyc, cycles and Q.");
	puts("  cyc=   : Number of cycles defined at one sigma (from - to +). Overwrites cycles & Q.");
	puts("  cycles=: Number of cycles defined at the half-power bandwidth. Overwrites Q.");
	puts("  Q=     : Ratio of central frequency to (half-power) bandwidth.");
	puts("  MexHat : Use the complex Mexican hat wavelet. Q = sqrt(1.5).");
	puts("           By default V = 2, b0 = 0.5, s0 = 1 and J = 6.");
	puts("           It ignores Q, cyc, cycles and w0 parameters.");
	puts("");
	puts("Other parameters (not set by default).");
	puts("  fold    : Sum positive and negative (lag-)times (used in seismic interferometry).");
	puts("  TwoStage: Two-stage ts-PWS: divide dataset in 10 groups, linear stack per group and ");
	puts("            ts-PWS of the 10 stacked traces. If TwoStage=X use X groups instead of 10.");
	puts("  unbiased: Compensate for PWS bias (1/N) whenever wu=2.");
	puts("  bin     : Provide the input traces in a single binary file (\"filename\").");
	puts("            Use sac2bin code to collect all sac files into a single binary file.");
	puts("  osac=   : Modify the default name of the output files by adding the text given.");
	puts("");
	puts("DEFAULTS:");
	puts("  ts-PWS with a power of 2 (wu=2) using a Morlet wavelet with: ");
	puts("     w0 = pi * sqrt(2/ln 2) (i.e, Q=pi/sqrt(2) or cycles=sqrt(2)),");
	puts("     V = 4, b0 = 1 and s0 = 2.");
	puts("  Two-stage off, unbiased off, rm off and fold off.");
	puts("  Output files:"); 
	puts("     Linear stack: tl.sac");
	puts("     ts-PWS:       ts_pws.sac"); 
	puts("");
	puts("AUTHOR: Sergi Ventosa, 17/01/2018");
	puts("");
	puts("EXAMPLES");
	puts("  1) Minimal: Stack of the traces listed in filelist.txt using the time-scale PWS method");
	puts("              using the defaults defined above.");
	puts("       ts_pws1e filelist.txt");
	puts("  2) Same as 1 but reducing the frequency range to 3 octaves starting at 4 mHz and");
	puts("     applying folding (summing positive and negative lag times).");
	puts("       ts_pws1e filelist.txt osac=\"example2\" wu=2 rm fold fmin=0.004 J=3 verbose");
	puts("  3) Same as 2 but providing the traces in a single binary file and using the two-stage");
	puts("     stacking method. sac2bin, stores the sac traces into a single binary file \"data.bin\".");
	puts("       ts_pws1e data.bin osac=\"twostage\" wu=2 rm bin fold TwoStage=10 unbiased fmin=0.004 J=3 verbose");
	puts("  Examples 2 & 3 reproduce the traces shown in Fig. 5c of Ventosa et al (GJI 2017).");
	puts(""); 
	puts("Please, do not hesitate to send bugs, comments or improvements to sergiventosa(at)hotmail.com\n");
}

/*****************************************************/
/* Functions to read the data and write the results. */ 
/*****************************************************/
int ReadData (float **sigall, time_t **time, t_hdr *hdr, char *filein, int bin, int verbose) {
	float *sig, beg, dt, beg1, dt1;
	unsigned int nsmp, nsamp, itr, mtr, nskip;
	int nerr=0, max, ia1;
	char **filenames = NULL, *filename=NULL;
	FILE *fid;
	size_t nitems;
	t_ccheader binhdr;
	
	if (filein == NULL) {printf("tspws_main: NULL filename\n"); return -1;}
	if (NULL == (fid = fopen(filein, "r"))) {
		printf("tspws_main: cannot open the %s file\n", filein);
		return -2;
	}
	
	/* Read the header. */
	if (bin) {
		nitems = fread (&binhdr, sizeof(t_ccheader), 1, fid);
		if (nitems != 1) 
			printf ("tspws_main: %s is shorter than predicted (header).\n", filein);
		hdr->max  = binhdr.nlags;
		hdr->mtr  = binhdr.nseq;
		hdr->evla = binhdr.stlat1;
		hdr->evlo = binhdr.stlon1;
		hdr->stla = binhdr.stlat2;
		hdr->stlo = binhdr.stlon2;
		hdr->dt   = (binhdr.lag2 - binhdr.lag1)/(float)(binhdr.nlags - 1);
		hdr->beg  = binhdr.lag1;
	} else {
		if ( (nerr = CreateFilelist (&filenames, &hdr->mtr, filein)) ) {
			printf("tspws_main: cannot read the %s file (CreateFileList error = %d)\n", filein, nerr);
			return -2;
		}
		if (hdr->mtr > 0) {
			/* Read the header of the first file. */
			filename = filenames[0];
			rsach (filename, &nerr, strlen(filename));            if (nerr) error_header (filename, nerr);
			sac_warning_off ();  /* Avoids lots of warnings when fields like stdp & knetwk are undefined. */
			getnhv ("npts",  &hdr->max,  &nerr, strlen("npts"));  if (nerr) return error_missing_info (filename, "npts", nerr); 
			getfhv ("stla",  &hdr->stla, &nerr, strlen("stla"));  if (nerr) warning_missing_info (filename, "stla", verbose);
			getfhv ("stlo",  &hdr->stlo, &nerr, strlen("stlo"));  if (nerr) warning_missing_info (filename, "stlo", verbose);
			getfhv ("evla",  &hdr->evla, &nerr, strlen("evla"));  if (nerr) warning_missing_info (filename, "evla", verbose);
			getfhv ("evlo",  &hdr->evlo, &nerr, strlen("evlo"));  if (nerr) warning_missing_info (filename, "evlo", verbose);
			getfhv ("delta", &hdr->dt,   &nerr, strlen("delta")); if (nerr) return error_missing_info (filename, "delta", nerr); 
			getfhv ("b",     &hdr->beg,  &nerr, strlen("b"));     if (nerr) return error_missing_info (filename, "b", nerr);
		} else {
			printf("tspws_main: nothing to do, %s is empty!\n", filein);
			return 0;
		}
	}
	max = hdr->max;
	mtr = hdr->mtr;
	beg = hdr->beg;
	dt  = hdr->dt;
	
	/* Allocated memory */
	if (NULL == (*sigall = (float *)mycalloc(mtr*max, sizeof(float)) )) nerr = 4;
	if (NULL == (*time   = (time_t *)mycalloc(mtr, sizeof(time_t)) ))   nerr = 4;
	if (nerr == 4) {
		myfree(*sigall);
		myfree(*time);
		printf ("tspws_main: Out of memory when reading %s (mtr = %d, npts = %d)\n", filein, mtr, max); 
		return nerr; 
	}
	
	nsamp = max;
	
	/* Reading and checking. */
	if (bin) {
		/* fseek (fid, mtr*sizeof(time_t), SEEK_CUR); */
		nitems = fread (*time, sizeof(time_t), mtr, fid);
		if (nitems != mtr) 
			printf ("tspws_main: %s is shorter than predicted (time).\n", filein);
		nitems = fread (*sigall, sizeof(float), mtr*max, fid);
		if (nitems != mtr*max) 
			printf ("tspws_main: %s is shorter than predicted (data).\n", filein);
		fclose (fid);
	} else {
		nsmp = nsamp;
		dt1  = dt;
		beg1 = beg;
		nskip = 0;
		for (itr=0; itr<mtr; itr++) {
			sig = *sigall + (itr-nskip)*max;
			filename = filenames[itr];
			rsac1(filename, sig, &ia1, &beg, &dt, &max, &nerr, strlen(filename));
			nsamp = (unsigned)ia1;
			
			if (nerr) {
				printf("tspws_main: Error reading %s file (rsac1, nerr=%d)\n", filename, nerr);
				return -2;
			}
			if (nsamp > nsmp)
				printf("tspws_main: WARNING: use only %u samples on %u trace\n", nsmp, itr);
			else if (nsamp < nsmp)
				printf("tspws_main: WARNING: trace %u has only %u samples\n", itr, nsamp);
			if (fabs(dt-dt1) > dt1*0.01) {
				printf("tspws_main: WARNING: trace %u has a different dt !\n", itr);
				printf("tspws_main: WARNING: skipping trace %u\n", itr);
				nskip += 1;
				continue;
			}
			if (fabs(beg1-beg) > dt1) {
				printf("tspws_main: WARNING: trace %u has a different beg !\n", itr);
				printf("tspws_main: WARNING: skipping trace %u\n", itr);
				nskip += 1;
				continue;
			}
		}
		mtr = itr - nskip;
		DestroyFilelist(filenames);
	}
	
	return nerr;
}

int wrsac(char *filename, float *y, t_hdr *hdr, char *kstnm, float user0) {
	float *dummy;
	/* float e, o; */
	int nerr;
	
	/* e = hdr->beg + (float)(hdr->nsamp-1)*hdr->dt;
	o = 0.; */
	
	if (NULL == (dummy = (float *)malloc(hdr->max*sizeof(float)) )) return -1;
	
	newhdr();
	setnhv ("npts",  &hdr->max,  &nerr, strlen("npts"));
	setfhv ("delta", &hdr->dt,   &nerr, strlen("delta"));
	setkhv ("kstnm",  kstnm,     &nerr, strlen("kstnm"), strlen(kstnm));
	setfhv ("stla",  &hdr->stla, &nerr, strlen("stla"));
	setfhv ("stlo",  &hdr->stlo, &nerr, strlen("stlo"));
	setfhv ("evla",  &hdr->evla, &nerr, strlen("evla"));
	setfhv ("evlo",  &hdr->evlo, &nerr, strlen("evlo"));
	setfhv ("user0", &user0,     &nerr, strlen("user0"));
	/* setfhv ("e",     &e,     &nerr, strlen("e")); */
	/* setfhv ("o",     &o,     &nerr, strlen("o")); */
	setfhv ("b",     &hdr->beg,  &nerr, strlen("b"));
	
	wsac0(filename, dummy, y, &nerr, strlen(filename));
	free(dummy);
	if (nerr) printf("\a wrsac: Error writing %s file\n", filename);
	return nerr;
}

void DestroyFilelist(char *p[]) {
	if (p) {
		if (p[0]) myfree(p[0]);
		myfree(p);
	}
}

int CreateFilelist (char **filename[], unsigned int *Tr, char *filelist) {
	long pos;
	unsigned int n, N, er=0;
	char *str0, *str1, **files, *mem_filenames;
	FILE *fid;
	
	/* Test parameters. */
	if (!filelist) return 1;
	
	/* Open the file and check its size. */
	if (NULL == (fid = fopen(filelist, "r"))) 
		{ printf("\a CreateFilelist: cannot open %s file\n", filelist); return 3; }
	
	/* Allocated memory to read file names. */
	if (NULL == (str0 = (char *)mymalloc(1024*sizeof(char)) )) 
		{ printf("\a CreateFilelist: out of memory."); return 2; }
	
	/* Get the number of lines & filesize */
	for (N=0; fgets(str0, 1024, fid); N++);
	fseek(fid, 0L, SEEK_END); /* Not actually needed. */
	pos=ftell(fid);
	fseek(fid, 0L, SEEK_SET);
	
	/* Allocate memory */
	files = (char **)mycalloc(N, sizeof(char *));
	mem_filenames = (char *)mymalloc((pos + N)*sizeof(char *));
	*filename = files;
	
	/* Read the file */
	if (files == NULL || mem_filenames == NULL) {
		if ( !N ) 
			printf("\a CreateFilelist: %s is empty!\n", filelist);
		else {
			er = 2;
			printf("\a CreateFilelist: out of memory.");
			myfree(mem_filenames);
			myfree(files);
		}
		*Tr = 0;
	} else {
		pos = 0;
		for (n=0; n<N; n++) {
			if (NULL == fgets(str0, 1024, fid)) break;
			str1 = strchr(str0,'\n');        /* Find newline character                 */
			if (str1) *str1 = '\0';          /* to properly set the end of the string. */
			
			mem_filenames[pos] = '\0';     /* Needed on the strcat below if no folder. */
			files[n] = &mem_filenames[pos];
			strcat(files[n], str0);        /* Add the filename.                        */
			pos += strlen(files[n]) + 1;   /* Update position at the filenames memory. */
		}
		*Tr = N;
	}
	fclose(fid);
	
	myfree(str0);
	return er;
}

void strcat3 (char *sout, const char *s1, const char *s2, const char *s3) {
	strcpy(sout, s1);
	strcat(sout, s2);
	strcat(sout, s3);
}
/*********************************************************/
/* Low-level functions used when reading the parameters. */
/*********************************************************/
int isuint (const char *str) {
	unsigned int n, n0, N;

	N = strlen(str);
	for (n=0; str[n]==' '; n++);
	if (str[n]=='+') n++;
	for (n0=n; isdigit(str[n]); n++);
	return (n==N && n>n0) ? 1 : 0;
}

int isdouble (const char *str) {
	unsigned int n, N;

	N = strlen(str);
	for (n=0;  str[n]==' '; n++);
	if (str[n]=='+' || str[n] =='-') n++;
	for (; isdigit(str[n]); n++);
	if (str[n]=='.') 
		for (n++; isdigit(str[n]); n++);
	return (n==N) ? 1 : 0;
}

int RDuint (unsigned int* const x, const char *str) {
	if (!isuint(str)) return 1;
	*x = (unsigned)atoi(str);
	return 0;
}

int RDdouble (double* const x, const char *str) {
	if (!isdouble(str)) return 1;
	*x = atof(str);
	return 0;
}
