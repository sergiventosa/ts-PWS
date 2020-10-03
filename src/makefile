CC=gcc
DEBUG=-g -DDEBUG -pg
OPT = -O3 -march=native
CFLAGS=-c -Wall $(OPT) -ITools -IFWTa -I$(SACHOME)/include -DNoThreads -fopenmp
LFLAGS=-Wall $(OPT) -fopenmp
CLIBS=-lm

libsac = $(SACHOME)/lib/sacio.a

VPATH = . ./Tools ./FWTa

all: ts_pws sac2bin

ts_pws: ts_pws1f.o ts_pws1f_lib.o wavelet_v7.o wavelet_def_v7.o wavelet_mem_v7.o cdotx.o myallocs.o prnmsg.o
	$(CC) $(LFLAGS) -o ts_pws ts_pws1f.o ts_pws1f_lib.o wavelet_v7.o wavelet_def_v7.o wavelet_mem_v7.o cdotx.o myallocs.o prnmsg.o $(libsac) $(CLIBS)
	
ts_pws1f_lib.o: ts_pws1f_lib.c ts_pws1f_lib.h
	$(CC) $(CFLAGS) ts_pws1f_lib.c
	
ts_pws1f.o: ts_pws1f.c ts_pws1f_lib.h
	$(CC) $(CFLAGS) ts_pws1f.c
	
sac2bin: sac2bin.o
	$(CC) $(LFLAGS) -o sac2bin sac2bin.o $(libsac) $(CLIBS)
	
install: ts_pws sac2bin
	mkdir -p ../bin
	install -s ts_pws sac2bin ../bin
	
clean:
	rm -rf *o ts_pws sac2bin

