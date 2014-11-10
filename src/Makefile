CC=g++
CFLAGS=-c -Wall
INCLUDES= -I../include -I../library/seqPatch/src/

all: pileup prior detect

pileup: qDNAmod_pileup.cpp 
	$(CC) $(INCLUDES) -o ../bin/qDNAmod_pileup qDNAmod_pileup.cpp

prior: qDNAmod_prior.cpp
	$(CC) $(INCLUDES) -o ../bin/qDNAmod_prior qDNAmod_prior.cpp

detect: qDNAmod_detect.cpp
	$(CC) $(INCLUDES) -o ../bin/qDNAmod_detect qDNAmod_detect.cpp
clean:
	rm -rf ../bin/*.o ../bin/qDNAmod_pileup ../bin/qDNAmod_prior ../bin/qDNAmod_detect


