CC=g++
CFLAGS=-c -Wall
INCLUDES= -I../include -I../library/seqPatch/src/

all: pileup

pileup: seqBMA_pileup.cpp 
	$(CC) $(INCLUDES) -o ../bin/seqBMA_pileup seqBMA_pileup.cpp
clean:
	rm -rf ../bin/*.o ../bin/seqBMA_pileup


