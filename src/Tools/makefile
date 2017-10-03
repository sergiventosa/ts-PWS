
all: myallocs.o prnmsg.o

myallocs.o: myallocs.c myallocs.h
	$(CC) $(CFLAGS) myallocs.c

prnmsg.o: prnmsg.c prnmsg.h
	$(CC) $(CFLAGS) prnmsg.c

clean:
	rm -rf *o

