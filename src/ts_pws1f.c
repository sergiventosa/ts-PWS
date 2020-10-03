/*****************************************************************************/
/* Program to perform time-scale phase-weighted stack (ts-PWS) using frames  */
/* of continuous wavelets, including the two-stage stack and the unbiased    */
/* phase coherence strategies. (Ventosa et al., GJI, 2017)                   */
/*                                                                           */
/*  Ventosa, S., Schimmel, M., & Stutzmann, E., 2017. Extracting surface     */
/*  waves, hum and normal modes: Time-scale phase-weighted stack and beyond, */
/*  Geophysical Journal International, 211(1), 30-44, doi:10.1093/gji/ggx284 */
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
/* **** 2019 ****                                                            */
/* Feb20 (1f) New feature:                                                   */
/*            - Subsampled linear stacks and ts-PWS.                         */
/*            - Added optimal weights per trace.                             */
/* Aug31 (1f) Bug corrections:                                               */
/*            - Set the lcalda flag. This way SAC computes dist, az, baz and */
/*            gcarc from station and event coordinates. Otherwise, this is   */
/*            done when the sac-file was read by SAC.                        */
/* Sep10 (1f) Much better openMP parallization of TwoStage_subsmpl_float().  */
/*            Bug corrections:                                               */
/*            - Solved memory error in tspws_subsmpl_float() that appeared   */
/*            when subsampling a very small number of sequences.             */
/* **** 2020 ****                                                            */
/* Sep16 (1f) Modified the format of the input binary file (msacs) to save   */
/*            the actual "zero" lag time for each correlation pair.          */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sacio.h>
#include <math.h>

#include "sac2bin.h"
#include "myallocs.h"
#include "ts_pws1f_lib.h"

void infooo();
void usage();

/* Caution: not checking the final result or intermediate results fits the range of unsigned int. */
unsigned int binomial_coefficient (unsigned int n, unsigned int d) {
	unsigned int i, out;
	
	if (d > n-d) d = n-d;
	out = 1;
	for (i=n-d+1; i<=n; i++) out *= i;
	for (i=2; i<=d; i++) out /= i;
	
	return out;
}

int ReadData (float **sigall, time_t **time, float **lag0, t_hdr *hdr, char *filein, int bin, int verbose);
int wrsac(char *filename, float *y, t_hdr *hdr, char *kinst, float user0);
int wrbin(char *filename, float **y, unsigned int M, unsigned int n, t_hdr *hdr, char *kinst, unsigned int *mtr);

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

void Clean_tsPWSout (t_tsPWS_out *out) {
	myfree(out->ls);
	myfree(out->tsPWS);
	myfree(out->ls_sim);
	myfree(out->tsPWS_sim);
	myfree(out->ls_misfit);
	myfree(out->tsPWS_misfit);
	myfree(out->ls_steps);
	myfree(out->tsPWS_steps);
	if (out->ls_subsmpl) myfree(out->ls_subsmpl[0]);
	myfree(out->ls_subsmpl);
	if (out->tsPWS_subsmpl) myfree(out->tsPWS_subsmpl[0]);
	myfree(out->tsPWS_subsmpl);
	memset(out, '\0', sizeof(*out));
}

/********************************************************************************************/
/* Main function: Reads the parameters and the data, allocates memory and save the results. */
/********************************************************************************************/
int main(int argc, char *argv[]) {
	t_tsPWS tspws = {-1, 0, 0, 4, 2., 1.0, PI*sqrt(2/log(2)), 2., 0., 0., 2., 
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
				NULL, NULL, NULL, NULL};  /* Default parameters. */
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
		else if (!strncmp(argv[i], "bin",          3)) tspws.bin      =  1; /* True */
		else if (!strncmp(argv[i], "verbose",      7)) tspws.verbose  =  1; /* True */
		else if (!strncmp(argv[i], "fold",         4)) tspws.fold     =  1; /* True */
		else if (!strncmp(argv[i], "unbiased",     8)) tspws.unbiased =  1; /* True */
		else if (!strncmp(argv[i], "MexHat",       6)) tspws.type     = -3;
		else if (!strncmp(argv[i], "AllSteps",     8)) tspws.AllSteps =  1; /* True */
		else if (!strncmp(argv[i], "subsmpl_N=",  10)) er = RDuint(&tspws.subsmpl_N, argv[i] + 10);
		else if (!strncmp(argv[i], "subsmpl_prob=",13)) er = RDdouble(&tspws.subsmpl_p, argv[i] + 13);
		else if (!strncmp(argv[i], "jackknife_n=", 12)) er = RDuint(&tspws.jackknife_n, argv[i] + 12);
		else if (!strncmp(argv[i], "jackknife_d=", 12)) er = RDuint(&tspws.jackknife_d, argv[i] + 12);
		else if (!strncmp(argv[i], "obin",          4)) tspws.obin = 1;
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
		} else if (!strncmp(argv[i], "osac=", 5)) {
			tspws.fileout = argv[i] + 5;
		 	if (!strlen(tspws.fileout)) tspws.fileout = NULL;
		} else if (!strncmp(argv[i], "kinst=", 6)) {
			tspws.kinst  = argv[i] + 6;
			tspws.lkinst = 1; /* True */
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
	er = ReadData (&in.sigall, &in.time, &in.lag0, &in.hdr, tspws.filein, tspws.bin, tspws.verbose);
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
	/* Memory for convergence mesurements */
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
	/* Memory for subsampling */
	/*   jackknife            */
	if (tspws.jackknife_n) {
		tspws.subsmpl_N = 0;
		tspws.subsmpl_p = 0;
		if (tspws.jackknife_d == 0 || tspws.jackknife_d >= tspws.jackknife_n) { tspws.jackknife_d = 0; tspws.jackknife_n = 0; }
		else out.M = binomial_coefficient (tspws.jackknife_n, tspws.jackknife_d);
	}
	/*   random subset selection  */
	if (tspws.subsmpl_N) {
		if (tspws.subsmpl_p < 0 || tspws.subsmpl_p > 1) { tspws.subsmpl_N = 0; tspws.subsmpl_p = 0; }
		else out.M = tspws.subsmpl_N;
	}
	
	if (tspws.jackknife_n || tspws.subsmpl_N) {
		if (NULL == (out.mtr_subsmpl = (unsigned int *)mycalloc(out.M, sizeof(unsigned int)) )) er = 4;
		if (NULL == (out.ls_subsmpl = (float **)mycalloc(out.M, sizeof(float *)) )) er = 4;
		if (!er) if (NULL == (out.ls_subsmpl[0] = (float *)mycalloc(out.M*in.hdr.max, sizeof(float)) )) er = 4;
		for (i=1; i<out.M; i++) out.ls_subsmpl[i] = out.ls_subsmpl[i-1] + in.hdr.max;
		
		if (NULL == (out.tsPWS_subsmpl = (float **)mycalloc(out.M, sizeof(float *)) )) er = 4;
		if (!er) if (NULL == (out.tsPWS_subsmpl[0] = (float *)mycalloc(out.M*in.hdr.max, sizeof(float)) )) er = 4;
		for (i=1; i<out.M; i++) out.tsPWS_subsmpl[i] = out.tsPWS_subsmpl[i-1] + in.hdr.max;
	}
	
	if (er == 4) {
		Clean_tsPWSout (&out);
		printf ("main: Out of memory when reading %s (npts = %d)\n", tspws.filein, in.hdr.max);
		return er; 
	}
	
	if (tspws.subsmpl_N)
		for (i=1; i<out.M; i++) out.mtr_subsmpl[i] = out.mtr * tspws.subsmpl_p;
	
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
		
		wrsac(filename, out.ls + n, &hdr, (tspws.lkinst) ? tspws.kinst : "t-lin", (float)out.mtr);
		
		if (tspws.verbose) {
			printf("Output files:\n");
			printf("  Linear stack: %s\n", filename);
		}
		
		/* TS-PWS */
		if (tspws.fileout) strcat3 (filename, "ts_pws_", tspws.fileout, ".sac");
		else strcpy(filename, "ts_pws.sac");
		
		wrsac(filename, out.tsPWS + n, &hdr, (tspws.lkinst) ? tspws.kinst : "ts_pws", (float)out.mtr);
		
		if (tspws.verbose) 
			printf("  ts-PWS:       %s\n", filename);
		
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
		/* Subsampled binary files. */
		if (tspws.jackknife_n || tspws.subsmpl_N) {
			char *subsmpl_str;
			unsigned int m;
			
			ua1 = 25;
			if (tspws.fileout) ua1 += strlen(tspws.fileout);
			subsmpl_str = (char *)mycalloc(ua1, sizeof(char));
			
			if (tspws.fileout) sprintf(subsmpl_str, "_%s_subsmpl", tspws.fileout);
			else strcpy(subsmpl_str, "_subsmpl");
			
			if (tspws.obin) {
				sprintf(filename, "tl%s.bin", subsmpl_str);
				wrbin(filename, out.ls_subsmpl, out.M, n, &hdr, (tspws.lkinst) ? tspws.kinst : "t-lin", out.mtr_subsmpl);
				sprintf(filename, "ts_pws%s.bin", subsmpl_str);
				wrbin(filename, out.tsPWS_subsmpl, out.M, n, &hdr, (tspws.lkinst) ? tspws.kinst : "ts_pws", out.mtr_subsmpl);
			} else {
				for (m=0; m<out.M; m++) {
					/* Subsampled linear stacks */
					sprintf(filename, "tl%s_%d.sac", subsmpl_str, m);
					wrsac(filename, out.ls_subsmpl[m] + n, &hdr, (tspws.lkinst) ? tspws.kinst : "t-lin", (float)out.mtr_subsmpl[m]);
					
					sprintf(filename, "ts_pws%s_%d.sac", subsmpl_str, m);
					wrsac(filename, out.tsPWS_subsmpl[m] + n, &hdr, (tspws.lkinst) ? tspws.kinst : "ts_pws", (float)out.mtr_subsmpl[m]);
				}
			}
			
			myfree(subsmpl_str);
		}
		
		myfree(filename);
	}
	/* clean */
	Clean_tsPWSout (&out);
	myfree(in.sigall);
	myfree(in.lag0);
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
	puts("Ventosa, S., Schimmel, M., & E. Stutzmann, 2017. Extracting surface waves, hum and normal modes: Time-scale phase-weighted stack and beyond, Geophysical Journal International, 211(1), 30-44, doi:10.1093/gji/ggx284"); 
	puts("AUTHOR: Sergi Ventosa Rahuet (sergiventosa(at)hotmail.com)");
	puts("Last modification: 01/10/2020\n");
}

void usage() {
	puts("\nCompute the time-scale phase-weighted stack (ts-PWS), including the two-stage stack and the");
	puts("unbiased phase coherence strategies to reduce signal attenuation and increase noise attenuation.");
	puts("\nUSAGE: ts_pws filename parameter_list");
	puts("  filename: One file name (SAC) per line, traces must have same begin (b),");
	puts("            number of samples (nsmpl) and sampling interval (dt).");
	puts("  parameters are optional and can be provided in arbitrary order without any");
	puts("  blanck around '='.");
	puts("");
	puts("Most commonly used parameters");
	puts("  cyc=   : Number of cycles defined at the one sigma (from -sigma to +sigma).");
	puts("           This is equivalent to w0=cyc*pi in the (default) Morlet wavelet.");
	puts("  fmin=  : Minimum frequency analyzed (Hz).");
	puts("  info   : write background and main references to screen.");
	puts("           Just type: ts_pws info.");
	puts("  kinst= : Optional string to be written on the instrument sac-header field (char*8).");
	puts("  rm     : remove the mean from input data.");
	puts("  verbose: Print the parameters used and files stacked.");
	puts("  wu=    : PWS power.");
	puts("");
	puts("Parameters controlling the discretization of the Continuous Wavelet Transform.");
	puts("  b0=    : Sampling at the scale 1.");
	puts("  J=     : Number of octaves.");
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
	puts("            This option can reduce signal attenuation and increase noise attenuation.");
	puts("  unbiased: Compensate for the bias in the phase coherence estimator (1/N) when wu=2.");
	puts("            This option can increase noise attenuation.");
	puts("  bin     : Provide the input traces in a single binary file (\"filename\").");
	puts("            Use the sac2bin code to collect all sac files into a single binary file or");
	puts("            the obin parameter in the fast phase cross-correlation code (FastPCC, v1.0.2),");
	puts("            https://github.com/sergiventosa/FastPCC.");
	puts("  osac=   : Modify the default name of the output files by adding the text given.");
	puts("");
	puts("Parameters to compute resampled stacks (not done by default).");
	puts("  Random resampling:");
	puts("    subsmpl_N=    : Number of subsampled stacks.");
	puts("    subsmpl_prob= : Probability of stacking each sequence.");
	puts("  Jackknife binning according to the julian day of the year (currently for the TwoStage only):");
	puts("    jackknife_n=  : Number of data bins.");
	puts("    jackknife_d=  : Number of deletions.");
	puts("  obin    : Returns resampled stacks in a single binary file binary file.");
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
	puts("AUTHOR: Sergi Ventosa, 01/10/2020");
	puts("Version 1.0.2");
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
int isstring (char *s, unsigned int N) {
	unsigned int n;
	
	for (n=0; n<N; n++) 
		if (!isalnum(s[n]) && !isspace(s[n]) && !ispunct(s[n]) && s[n] != 0 ) return n+1;
	return 0;
}

int check_header (t_hdr *h) {
	int n, out = 0;
	
	if (fabs(h->stla) > 90 || fabs(h->stlo) > 180 || fabs(h->evla) > 90 || fabs(h->evlo) > 180) out = 1;
	if (isstring(h->net1, 8) || isstring(h->sta1, 8) || isstring(h->chn1, 8) || isstring(h->loc1, 8) ) out = 1; 
	if (isstring(h->net2, 8) || isstring(h->sta2, 8) || isstring(h->chn2, 8) || isstring(h->loc2, 8) ) out = 1; 
	
	if (out) {
		printf("\n lat=%f\n lon=%f\n", h->stla, h->stlo);
		printf(" lat=%f\n lon=%f\n", h->evla, h->evlo);
		printf(" net1=%.8s\n sta1=%.8s\n loc1=%.8s\n chn1=%.8s\n", h->net1, h->sta1, h->loc1, h->chn1);
		printf(" net2=%.8s\n sta2=%.8s\n loc2=%.8s\n chn2=%.8s\n", h->net2, h->sta2, h->loc2, h->chn2);
		if (fabs(h->stla) > 90)  printf("stla is not fine\n");
		if (fabs(h->stlo) > 180) printf("stlo is not fine\n");
		if (fabs(h->evla) > 90)  printf("stla is not fine\n");
		if (fabs(h->evlo) > 180) printf("stlo is not fine\n");
		if ((n = isstring(h->net1, 8))) { n--; printf("net1[%d] is not fine (%d)\n", n, h->net1[n]); }
		if ((n = isstring(h->sta1, 8))) { n--; printf("sta1[%d] is not fine (%d)\n", n, h->sta1[n]); }
		if ((n = isstring(h->chn1, 8))) { n--; printf("chn1[%d] is not fine (%d)\n", n, h->chn1[n]); }
		if ((n = isstring(h->loc1, 8))) { n--; printf("loc1[%d] is not fine (%d)\n", n, h->loc1[n]); }
		if ((n = isstring(h->net2, 8))) { n--; printf("net2[%d] is not fine (%d)\n", n, h->net2[n]); }
		if ((n = isstring(h->sta2, 8))) { n--; printf("sta2[%d] is not fine (%d)\n", n, h->sta2[n]); }
		if ((n = isstring(h->chn2, 8))) { n--; printf("chn2[%d] is not fine (%d)\n", n, h->chn2[n]); }
		if ((n = isstring(h->loc2, 8))) { n--; printf("loc2[%d] is not fine (%d)\n", n, h->loc2[n]); }
	}
	return out;
}

int ReadData (float **sigall, time_t **time, float **lag0, t_hdr *hdr, char *filein, int bin, int verbose) {
	float *sig, beg, dt, beg1, dt1;
	unsigned int nsmp, nsamp, itr, mtr, nskip, n;
	int nerr=0, max, ia1;
	char **filenames = NULL, *filename=NULL;
	FILE *fid;
	size_t nitems;
	t_ccheader binhdr;
	
	if (filein == NULL) {printf("tspws_main: NULL filename\n"); return -1; }
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
		hdr->stel = binhdr.stel2;
		hdr->dt   = (binhdr.lag2 - binhdr.lag1)/(float)(binhdr.nlags - 1);
		hdr->beg  = binhdr.lag1;
		strncpy(hdr->net1, binhdr.net1, 9);
		strncpy(hdr->sta1, binhdr.sta1, 9);
		strncpy(hdr->loc1, binhdr.loc1, 9);
		strncpy(hdr->chn1, binhdr.chn1, 9);
		strncpy(hdr->net2, binhdr.net2, 9);
		strncpy(hdr->sta2, binhdr.sta2, 9);
		strncpy(hdr->loc2, binhdr.loc2, 9);
		strncpy(hdr->chn2, binhdr.chn2, 9);
		
		if (check_header (hdr) ) printf("\a ReadData: Found the header corrupted when reading the %s file\n", filein);
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
			getnhv ("npts",   &hdr->max,  &nerr, strlen("npts"));  if (nerr) return error_missing_info (filename, "npts", nerr);
			getfhv ("stla",   &hdr->stla, &nerr, strlen("stla"));  if (nerr) warning_missing_info (filename, "stla", verbose);
			getfhv ("stlo",   &hdr->stlo, &nerr, strlen("stlo"));  if (nerr) warning_missing_info (filename, "stlo", verbose);
			getfhv ("evla",   &hdr->evla, &nerr, strlen("evla"));  if (nerr) warning_missing_info (filename, "evla", verbose);
			getfhv ("evlo",   &hdr->evlo, &nerr, strlen("evlo"));  if (nerr) warning_missing_info (filename, "evlo", verbose);
			getfhv ("delta",  &hdr->dt,   &nerr, strlen("delta")); if (nerr) return error_missing_info (filename, "delta", nerr);
			getfhv ("b",      &hdr->beg,  &nerr, strlen("b"));     if (nerr) return error_missing_info (filename, "b", nerr);
			
			getfhv ("stel",   &hdr->stel, &nerr, strlen("stel"));
			getkhv ("knetwk",  hdr->net2, &nerr, strlen("knetwk"), strlen(hdr->net2));
			getkhv ("kstnm",   hdr->sta2, &nerr, strlen("kstnm"),  strlen(hdr->sta2));
			getkhv ("khole",   hdr->loc2, &nerr, strlen("khole"),  strlen(hdr->loc2));
			getkhv ("kcmpnm",  hdr->chn2, &nerr, strlen("kcmpnm"), strlen(hdr->chn2)); if (nerr) warning_missing_info (filename, "kcmpnm", verbose);
			
			getkhv ("kuser0",  hdr->net1, &nerr, strlen("kuser0"), strlen(hdr->net1));
			getkhv ("kevnm",   hdr->sta1, &nerr, strlen("kevnm"),  strlen(hdr->sta1));
			getkhv ("kuser1",  hdr->loc1, &nerr, strlen("kuser1"), strlen(hdr->loc1));
			getkhv ("kuser2",  hdr->chn1, &nerr, strlen("kuser2"), strlen(hdr->chn1)); if (nerr) warning_missing_info (filename, "kuser2", verbose);
			
			if (check_header (hdr) ) printf("\a ReadData: Found the header corrupted when reading the %s file\n", filename);
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
	if (NULL == (*lag0   = (float *)mycalloc(mtr, sizeof(float)) ))     nerr = 4;
	if (nerr == 4) {
		myfree(*sigall);
		myfree(*lag0);
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
		nitems = fread (*lag0, sizeof(float), mtr, fid);
		if (nitems != mtr) 
			printf ("tspws_main: %s is shorter than predicted (lag0).\n", filein);
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
				printf("tspws_main: WARNING: using only %u samples on %u trace\n", nsmp, itr);
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
	
	/* Check for zero traces */
	for (itr=0; itr<mtr; itr++) {
		sig = *sigall + itr*max;
		for (n=0; n<max; n++) if ( sig[n] != 0.) break;
		if (n == max) {
			printf("tspws_main: %s, trace %d of %d is ZERO\n", filein, itr, mtr);
			continue;
		}
	}
	
	return nerr;
}

int wrsac(char *filename, float *y, t_hdr *hdr, char *kinst, float user0) {
	float *dummy;
	int nerr, lcalda = 1;
	
	if (NULL == (dummy = (float *)malloc(hdr->max*sizeof(float)) )) return -1;
	
	if (check_header (hdr) ) printf("\a wrsac: Found the header corrupted when writing the %s file\n", filename);
	
	newhdr();
	setihv ("iztype", "io",  &nerr, strlen("iztype"), strlen("io"));
	setnhv ("npts",  &hdr->max,  &nerr, strlen("npts"));
	setfhv ("delta", &hdr->dt,   &nerr, strlen("delta"));
	setkhv ("kinst", kinst,      &nerr, strlen("kinst"), strlen(kinst));
	setfhv ("stla",  &hdr->stla, &nerr, strlen("stla"));
	setfhv ("stlo",  &hdr->stlo, &nerr, strlen("stlo"));
	setfhv ("stel",  &hdr->stel, &nerr, strlen("stel"));
	setfhv ("evla",  &hdr->evla, &nerr, strlen("evla"));
	setfhv ("evlo",  &hdr->evlo, &nerr, strlen("evlo"));
	setfhv ("user0", &user0,     &nerr, strlen("user0"));
	setfhv ("b",     &hdr->beg,  &nerr, strlen("b"));
	
	setkhv ("knetwk", hdr->net2, &nerr, strlen("knetwk"), (strlen(hdr->net2) < 8) ? strlen(hdr->net2) : 8 );
	setkhv ("kstnm",  hdr->sta2, &nerr, strlen("kstnm"),  (strlen(hdr->sta2) < 8) ? strlen(hdr->sta2) : 8 );
	setkhv ("khole",  hdr->loc2, &nerr, strlen("khole"),  (strlen(hdr->loc2) < 8) ? strlen(hdr->loc2) : 8 );
	setkhv ("kcmpnm", hdr->chn2, &nerr, strlen("kcmpnm"), (strlen(hdr->chn2) < 8) ? strlen(hdr->chn2) : 8 );
	
	setkhv ("kuser0", hdr->net1, &nerr, strlen("kuser0"), (strlen(hdr->net1) < 8) ? strlen(hdr->net1) : 8 );
	setkhv ("kevnm",  hdr->sta1, &nerr, strlen("kevnm"),  (strlen(hdr->sta1) < 8) ? strlen(hdr->sta1) : 8 );
	setkhv ("kuser1", hdr->loc1, &nerr, strlen("kuser1"), (strlen(hdr->loc1) < 8) ? strlen(hdr->loc1) : 8 );
	setkhv ("kuser2", hdr->chn1, &nerr, strlen("kuser2"), (strlen(hdr->chn1) < 8) ? strlen(hdr->chn1) : 8 );
	
	setlhv ("lcalda", &lcalda, &nerr, strlen("lcalda"));
	
	wsac0(filename, dummy, y, &nerr, strlen(filename));
	free(dummy);
	
	if (nerr) printf("\a wrsac: Error writing the %s file\n", filename);
	return nerr;
}

int wrbin (char *filename, float **y, unsigned int M, unsigned int n, t_hdr *hdr, char *kinst, unsigned int *mtr) {
	float *data;
	time_t *time;
	t_ccheader binhdr;
	FILE *fid;
	int nerr=0, NLAGS;
	unsigned int m;
	
	if (check_header (hdr) ) printf("\a wrsac: Found the header corrupted when writing the %s file\n", filename);
	
	if (filename == NULL) { printf("wrbin: NULL filename\n"); return -1; }
	
	/* Header */
	binhdr.nlags   = hdr->max;
	binhdr.nseq    = M;
	binhdr.stlat1  = hdr->evla;
	binhdr.stlon1  = hdr->evlo;
	binhdr.stlat2  = hdr->stla;
	binhdr.stlon2  = hdr->stlo;
	binhdr.stel2   = hdr->stel;
	binhdr.lag1    = hdr->beg;
	binhdr.tlength = hdr->dt * (float)(binhdr.nlags - 1);
	binhdr.lag2    = binhdr.lag1 + binhdr.tlength;
	strncpy(binhdr.method, kinst, 8);
	strncpy(binhdr.net1, hdr->net1, 8);
	strncpy(binhdr.sta1, hdr->sta1, 8);
	strncpy(binhdr.loc1, hdr->loc1, 8);
	strncpy(binhdr.chn1, hdr->chn1, 8);
	strncpy(binhdr.net2, hdr->net2, 8);
	strncpy(binhdr.sta2, hdr->sta2, 8);
	strncpy(binhdr.loc2, hdr->loc2, 8);
	strncpy(binhdr.chn2, hdr->chn2, 8);
	
	/* Allocated memory */
	NLAGS = hdr->max;
	if (NULL == (data = (float *)mycalloc(M*NLAGS, sizeof(float)) )) nerr = 4;
	if (NULL == (time = (time_t *)mycalloc(M, sizeof(time_t)) )) nerr = 4;
	if (nerr == 4) {
		myfree(data);
		myfree(time);
		printf ("wrbin: Out of memory when writing %s (mtr = %d, npts = %d)\n", filename, M, NLAGS); 
		return nerr; 
	}
	
	for (m=0; m<M; m++) time[m] = (time_t)mtr[m];
	for (m=0; m<M; m++) memcpy(data + m*NLAGS, y[m] + n, NLAGS*sizeof(float));
	
	if (NULL == (fid = fopen(filename, "w"))) {
		printf("tspws_main: cannot open the %s file\n", filename);
		return -2;
	}
	
	fwrite (&binhdr, sizeof(t_ccheader), 1, fid);
	fwrite (time, sizeof(time_t), M, fid);
	fwrite (data, sizeof(float), M*NLAGS, fid);
	
	fclose(fid);
	myfree(data);
	myfree(time);
	
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
