CC=g++
CFLAGS=-c -Wall
INCLUDES= -I../include -I../library/seqPatch/src/

all: pileup prior detect

pileup: seqBMA_pileup.cpp 
	$(CC) $(INCLUDES) -o ../bin/seqBMA_pileup seqBMA_pileup.cpp

prior: seqBMA_prior.cpp
	$(CC) $(INCLUDES) -o ../bin/seqBMA_prior seqBMA_prior.cpp

detect: seqBMA_detect.cpp
	$(CC) $(INCLUDES) -o ../bin/seqBMA_detect seqBMA_detect.cpp
clean:
	rm -rf ../bin/*.o ../bin/seqBMA_pileup ../bin/seqBMA_prior ../bin/seqBMA_detect


