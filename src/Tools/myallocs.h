#ifndef MYALLOCS_H
#define MYALLOCS_H

#ifndef NoThreads
int myallocs_init(void);
void myallocs_destroy(void);
#endif
void *mymalloc (int n);
void *mycalloc (int n, int f);
void *myrealloc (void *, int n);
void myfree (void *p);

#endif
