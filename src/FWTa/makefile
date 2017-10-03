VPATH = ../Tools

all: wavelet_v7.o wavelet_def_v7.o wavelet_mem_v7.o cdotx.o myallocs.o myarrays.o prnmsg.o

wavelet_v7.o: wavelet_v7.c wavelet_v7.h
	$(CC) $(CFLAGS) wavelet_v7.c

wavelet_def_v7.o: wavelet_def_v7.c wavelet_v7.h
	$(CC) $(CFLAGS) wavelet_def_v7.c

wavelet_mem_v7.o: wavelet_mem_v7.c wavelet_v7.h
	$(CC) $(CFLAGS) wavelet_mem_v7.c
	
cdotx.o: cdotx.c cdotx.h
	$(CC) $(CFLAGS) cdotx.c

clean:
	rm -rf *.o *.a

