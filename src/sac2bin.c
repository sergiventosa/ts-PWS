#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sacio.h>
#include "sac2bin.h"

int sac2bin_main (char *outfile, char *infiles);
void usage ();
time_t utc_mktime (struct tm *tm);
void DestroyFilelist(char *p[]);
int CreateFilelist (char **filename[], unsigned int *Tr, char *filelist);

int main(int argc, char *argv[]) {
	char *infiles, *outfiles;
	
	if (argc < 2) usage();
	else {
		infiles  = argv[1];
		outfiles = argv[2];
		sac2bin_main (outfiles, infiles);
	}
	
	return 0;
}

int nerr_print (char *filename, int nerr) {
	printf ("\a sac2bin_main: Error reading %s header (nerr = %d)\n", filename, nerr); 
	return 2;
}

void warning_missing_info (char *filename, char *nick) {
	printf ("\a sac2bin_main: %s is not defined in %s.\n", filename, nick); 
}

void warning_diff (char *filename, char *nick) {
	printf ("\a sac2bin_main: %s is not defined in %s or different than the one defined in the first sequence.\n", filename, nick); 
}


int sac2bin_main (char *outfile, char *infiles) {
	t_ccheader *hdr;
	float *data=NULL, *trace=NULL, *lag0=NULL;
	float stlat1, stlon1, stlat2, stlon2;
	float LAG1, LAG2, DT;
	float lag1, lag2, dt;
	unsigned int itr, mtr, nskip;
	int nerr, ia1, nlags, NLAGS, skiptime;
	int year, yday, hour, min, sec, msec;
	char **filenames = NULL, *filename=NULL;
	FILE *fid;
	time_t *time=NULL;
	struct tm tm;
	
	if (outfile == NULL || infiles == NULL) { 
		printf("\a sac2bin_main: NULL inputs\n"); 
		return -1; 
	}
	if (NULL == (fid = fopen(infiles, "r"))) {
		printf("\a sac2bin_main: cannot open %s file\n", infiles);
		return -2;
	}
	
	if ( (nerr = CreateFilelist (&filenames, &mtr, infiles)) ) {
		printf("\a sac2bin_main: cannot read %s file (CreateFileList error = %d)\n", infiles, nerr);
		return -2;
	}
	
	/* Read the header of the first file. */
	if (NULL == (hdr = (t_ccheader *)calloc(1, sizeof(t_ccheader)) )) {
		printf("\a sac2bin_main: Out of memory (hdr).\n");
		return -4;
	}
	
	filename = filenames[0];
	rsach (filename, &nerr, strlen(filename));       if (nerr) return nerr_print (filename, nerr);
	// Not supported from v102.0
	// sac_warning_off ();  /* Avoids lots of warnings when fields such as stdp or knetwk are undefined. */
	getkhv ("kinst",  hdr->method, &nerr, strlen("kinst"),  8); if (nerr) strcpy(hdr->method, "");
	getkhv ("kuser0", hdr->net1,   &nerr, strlen("kuser0"), 8); if (nerr) strcpy(hdr->net1, "");
	getkhv ("kevnm",  hdr->sta1,   &nerr, strlen("kevnm"),  8); if (nerr) strcpy(hdr->sta1, "");
	getkhv ("kuser1", hdr->loc1,   &nerr, strlen("kuser1"), 8); if (nerr) strcpy(hdr->loc1, "");
	getkhv ("kuser2", hdr->chn1,   &nerr, strlen("kuser2"), 8); if (nerr) strcpy(hdr->chn1, "");
	getkhv ("knetwk", hdr->net2,   &nerr, strlen("knetwk"), 8); if (nerr) strcpy(hdr->net2, "");
	getkhv ("kstnm",  hdr->sta2,   &nerr, strlen("kstnm"),  8); if (nerr) strcpy(hdr->sta2, "");
	getkhv ("khole",  hdr->loc2,   &nerr, strlen("khole"),  8); if (nerr) strcpy(hdr->loc2, "");
	getkhv ("kcmpnm", hdr->chn2,   &nerr, strlen("kcmpnm"), 8); if (nerr) strcpy(hdr->chn2, "");
	getfhv ("evla",  &hdr->stlat1, &nerr, strlen("evla"));  if (nerr) warning_missing_info (filename, "evla");
	getfhv ("evlo",  &hdr->stlon1, &nerr, strlen("evlo"));  if (nerr) warning_missing_info (filename, "evlo");
	getfhv ("stla",  &hdr->stlat2, &nerr, strlen("stla"));  if (nerr) warning_missing_info (filename, "stla");
	getfhv ("stlo",  &hdr->stlon2, &nerr, strlen("stlo"));  if (nerr) warning_missing_info (filename, "stlo");
	getfhv ("delta", &DT,          &nerr, strlen("delta")); if (nerr) return nerr_print (filename, nerr);
	getnhv ("npts",  &ia1,         &nerr, strlen("npts"));  if (nerr) return nerr_print (filename, nerr);
	getfhv ("b",     &hdr->lag1,   &nerr, strlen("b"));     if (nerr) return nerr_print (filename, nerr);
	hdr->nlags = (uint32_t)ia1;
	hdr->lag2  = hdr->lag1 + DT*(ia1-1);
	hdr->tlength = fabs(hdr->lag2 - hdr->lag1);
	
	/* Allocated memory */
	NLAGS = hdr->nlags;
	LAG1  = hdr->lag1;
	LAG2  = hdr->lag2;

	if (NULL == (time = (time_t *)calloc(mtr, sizeof(time_t)) )) nerr = 4;
	if (NULL == (lag0 = (float *)calloc(mtr, sizeof(float)) )) nerr = 4;
	if (NULL == (data = (float *)calloc(mtr*NLAGS, sizeof(float)) )) nerr = 4;
	if (nerr == 4) {
		free(hdr);
		free(time);
		free(lag0);
		free(data);
		printf ("\a sac2bin_main: Out of memory when reading %s (npts = %d)\n", filename, mtr*NLAGS); 
		return 4; 
	}
	
	/*************************/
	/* Reading and checking. */
	/*************************/
	nskip = 0;
	for (itr=0; itr<mtr; itr++) {
		skiptime = 0;
		trace = data + (itr-nskip)*NLAGS;
		filename = filenames[itr];
		rsac1(filename, trace, &nlags, &lag1, &dt, &NLAGS, &nerr, strlen(filename));
		// Not supported from v102.0
		// sac_warning_off ();
		lag2 = lag1 + DT*(nlags-1);
		if (nerr) return nerr_print (filename, nerr);
		getnhv ("nzyear", &year,   &nerr, strlen("nzyear")); if (nerr) skiptime = 1;
		getnhv ("nzjday", &yday,   &nerr, strlen("nzjday")); if (nerr) skiptime = 1;
		getnhv ("nzhour", &hour,   &nerr, strlen("nzhour")); if (nerr) skiptime = 1;
		getnhv ("nzmin",  &min,    &nerr, strlen("nzmin"));  if (nerr) skiptime = 1;
		getnhv ("nzsec",  &sec,    &nerr, strlen("nzsec"));  if (nerr) skiptime = 1;
		getnhv ("nzmsec", &msec,   &nerr, strlen("nzmsec")); if (nerr) skiptime = 1; 
		getfhv ("o",  &lag0[itr-nskip], &nerr, strlen("o")); if (nerr) lag0[itr-nskip] = 0;
		getfhv ("evla",   &stlat1, &nerr, strlen("evla"));
		if (nerr || stlat1 != hdr->stlat1) warning_diff (filename, "evla");
		getfhv ("evlo",   &stlon1, &nerr, strlen("evlo"));
		if (nerr || stlon1 != hdr->stlon1) warning_diff (filename, "evlo");
		getfhv ("stla",   &stlat2, &nerr, strlen("stla"));
		if (nerr || stlat2 != hdr->stlat2) warning_diff (filename, "stla");
		getfhv ("stlo",   &stlon2, &nerr, strlen("stlo"));
		if (nerr || stlon2 != hdr->stlon2) warning_diff (filename, "stlo");
		
		if (!skiptime) {
			tm.tm_year = year-1900;
			tm.tm_mon  = 0;
			tm.tm_mday = yday;
			tm.tm_hour = hour;
			tm.tm_min  = min;
			tm.tm_sec  = sec;
			time[itr-nskip] = utc_mktime(&tm);
		}
		
		if (NLAGS < nlags)
			printf("\a sac2bin_main: WARNING: use only %u samples on %u trace\n", NLAGS, itr);
		else if (NLAGS > nlags)
			printf("\a sac2bin_main: WARNING: trace %u has only %u samples\n", itr, nlags);
		
		if (fabs(dt-DT) > DT*0.01) {
			printf("\a sac2bin_main: WARNING: trace %u has a different lag rate !\n", itr);
			nerr = 1;
		}
		if (fabs(LAG1-lag1) > DT) {
			printf("\a sac2bin_main: WARNING: trace %u has a different lag1 !\n", itr);
			nerr = 1;
		}
		if (fabs(LAG2-lag2) > DT) {
			printf("\a sac2bin_main: WARNING: trace %u has a different lag2 !\n", itr);
			nerr = 1;
		}
		if (nerr) {
			printf("\a sac2bin_main: WARNING: skipping trace %u\n", itr);
			nskip += 1;
			nerr = 0;
			continue;
		}
	}
	mtr = itr - nskip;
	hdr->nseq = mtr; 
	DestroyFilelist(filenames);

	/*******************************/
	/* Writing to the binary file .*/
	/*******************************/
	if (NULL == (fid = fopen(outfile, "w"))) {
		printf("\a sac2bin_main: cannot create %s file\n", outfile);
		return -2;
	}
	
	fwrite (hdr, sizeof(t_ccheader), 1, fid);
	fwrite (time, sizeof(time_t), mtr, fid);
	fwrite (lag0, sizeof(float), mtr, fid);
	fwrite (data, sizeof(float), mtr*NLAGS, fid);
	fclose (fid);
	
	free(hdr);
	free(time);
	free(lag0);
	free(data);
	
	return 0;
}

void usage () {
	puts("\nUSAGE: sac2bin \"List of sac files\" \"Output file name\"");
}

char *set_utc () {
	char *tz;
	
	tz = getenv ("TZ");
	if (tz) tz = strdup (tz);
	setenv ("TZ", "", 1);
	tzset ();

	return tz;
}

void restor_tz (char *tz) {
	if (tz) {
		setenv ("TZ", tz, 1);
		free (tz);
	} else unsetenv ("TZ");
	tzset ();
	free (tz);
}

time_t utc_mktime (struct tm *tm) {
	char *tz;
	time_t t;
	
	tz = set_utc ();
	t = mktime (tm);
	restor_tz (tz);
	
	return t;
}

void DestroyFilelist(char *p[]) {
	if (p) {
		if (p[0]) free(p[0]);
		free(p);
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
	if (NULL == (str0 = (char *)malloc(1024*sizeof(char)) )) 
		{ printf("\a CreateFilelist: out of memory."); return 2; }
	
	/* Get the number of lines & filesize */
	for (N=0; fgets(str0, 1024, fid); N++);
	fseek(fid, 0L, SEEK_END); /* Not actually needed. */
	pos=ftell(fid);
	fseek(fid, 0L, SEEK_SET);
	
	/* Allocate memory */
	files = (char **)calloc(N, sizeof(char *));
	mem_filenames = (char *)malloc((pos + N)*sizeof(char *));
	*filename = files;
	
	/* Read the file */
	if (files == NULL || mem_filenames == NULL) {
		er = 2;
		printf("\a CreateFilelist: out of memory.");
		free(mem_filenames);
		free(files);
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
	
	free(str0);
	return er;
}
