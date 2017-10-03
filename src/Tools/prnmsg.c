#include <stdlib.h>
#include <stdarg.h>

#ifndef MATLAB
#include <stdio.h>
#else 
#include "mex.h"
#endif

#ifdef NoThreads
/***********************/
/* No threaded version */
/***********************/

int prerror(const char *format, ...) {
	va_list args;
	char str[200];

	va_start(args, format);
	vsnprintf(str, 200, format, args);
	va_end(args);

#ifdef MATLAB
	mexPrintf("Error: %s\n",str);
#else
	printf("Error: %s\n",str);
#endif

	return 0;
}

#else /* NoThreads */

#ifdef _MSC_VER
/*************************/
/* Windows (MSC) version */
/*************************/

#include <windows.h>

HANDLE Busy;

void prerror_init (void) {
	Busy = CreateMutex(NULL, 0, NULL);
}

void prerror_destroy (void) {
	CloseHandle(Busy);
}

int prerror(const char *format, ...) {
	va_list args;
	char str[200];

	va_start(args, format);
	vsnprintf(str, 200, format, args);
	va_end(args);

	WaitForSingleObject(Busy, INFINITE);

#ifdef MATLAB
	mexPrintf("Error: %s\n",str);
#else
	printf("Error: %s\n",str);
#endif

	ReleaseMutex(Busy);

	return 0;
}

#else /* _MSC_VER */
/****************/
/* GNUC version */
/****************/

#include <semaphore.h>

sem_t busy;

void prerror_init (void) {
	sem_init(&busy, 0, 1);
}

void prerror_destroy (void) {
	sem_destroy(&busy);
}

int prerror(const char *format, ...) {
	va_list args;
	char str[200];

	va_start(args, format);
	vsnprintf(str, 200, format, args);
	va_end(args);

	sem_wait(&busy);

#ifdef MATLAB
	mexPrintf("Error: %s\n",str);
#else
	printf("Error: %s\n",str);
#endif

	sem_post(&busy);

	return 0;
}

#endif /* _MSC_VER */

#endif /* NoThreads */
