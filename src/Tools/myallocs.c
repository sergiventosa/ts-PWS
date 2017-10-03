#include <stdlib.h>

#ifdef MATLAB
#include "matrix.h"
#endif

#ifdef NoThreads
#include <string.h>
/***********************/
/* No threaded version */
/***********************/

void *mymalloc(int n) {
	void *p;

#ifdef MATLAB
	p = mxMalloc(n); 
#else
	p = malloc(n);
#endif

	return p;
}

void *myrealloc(void *p, int n) {
#ifdef MATLAB
	p = mxRealloc(p, n); 
#else
	p = realloc(p, n);
#endif

	return p;
}

void *mycalloc(int n, int f) { 
	void *p;

#ifdef MATLAB
	p = mxCalloc(n, f);
#else
	p = calloc(n, f);
#endif

	return p;
}

void myfree(void *p) { 
#ifdef MATLAB
	mxFree(p);
#else
	free(p);
#endif
}

#else /* NoThreads */

#ifdef _MSC_VER
/*************************/
/* Windows (MSC) version */
/*************************/

#include <windows.h>

HANDLE NoMemOp;

int myallocs_init (void) {
	if (NULL == (NoMemOp = CreateMutex(NULL, 0, NULL) )) return -1;
	else return 0;
}

void myallocs_destroy (void) {
	CloseHandle(NoMemOp);
}

void *mymalloc(int n) {
	void *p;

	WaitForSingleObject(NoMemOp, INFINITE);
#ifdef MATLAB
	p = mxMalloc(n); 
#else
	p = malloc(n);
#endif
	ReleaseMutex(NoMemOp);

	return p;
}

void *myrealloc(void *p, int n) {
	WaitForSingleObject(NoMemOp, INFINITE);
#ifdef MATLAB
	p = mxRealloc(p, n); 
#else
	p = realloc(p, n);
#endif
	ReleaseMutex(NoMemOp);

	return p;
}

void *mycalloc(int n, int f) { 
	void *p;

	WaitForSingleObject(NoMemOp, INFINITE);
#ifdef MATLAB
	p = mxCalloc(n, f);
#else
	p = calloc(n, f);
#endif
	ReleaseMutex(NoMemOp);

	return p;
}

void myfree(void *p) { 
	WaitForSingleObject(NoMemOp, INFINITE);
#ifdef MATLAB
	mxFree(p);
#else
	free(p);
#endif
	ReleaseMutex(NoMemOp);
}

#else /* _MSC_VER */
/****************/
/* GNUC version */
/****************/

#include <semaphore.h>

sem_t NoMemOp;

int myallocs_init (void) {
	sem_init(&NoMemOp, 0, 1);
	return 0;
}

void myallocs_destroy (void) {
	sem_destroy(&NoMemOp);
}

void *mymalloc(int n) {
	void *p;

	sem_wait(&NoMemOp);
#ifdef MATLAB
	p = mxMalloc(n); 
#else
	p = malloc(n);
#endif
	sem_post(&NoMemOp); 

	return p;
}

void *myrealloc(void *p, int n) {
	sem_wait(&NoMemOp);
#ifdef MATLAB
	p = mxRealloc(p, n); 
#else
	p = realloc(p, n);
#endif
	sem_post(&NoMemOp); 

	return p;
}

void *mycalloc(int n, int f) { 
	void *p;

	sem_wait(&NoMemOp);
#ifdef MATLAB
	p = mxCalloc(n, f);
#else
	p = calloc(n, f);
#endif
	sem_post(&NoMemOp);

	return p;
}

void myfree(void *p) { 
	sem_wait(&NoMemOp);
#ifdef MATLAB
	mxFree(p);
#else
	free(p);
#endif
	sem_post(&NoMemOp);
}

#endif /* _MSC_VER */

#endif /* NoThreads */

