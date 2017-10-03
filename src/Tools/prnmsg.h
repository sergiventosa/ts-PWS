#ifndef PRNMSG_H
#define PRNMSG_H

#ifndef NoThreads
void prerror_init(void);
void prerror_destroy(void);
#endif
int prerror(const char *args, ...);

#endif
