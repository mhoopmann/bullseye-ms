#Set these variables if needed
C = gcc
CC = g++
FLAGS = -O3 -static -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILEOFFSET_BITS=64 -DGCC
MSTOOLKIT = ../MSToolkit
INCLUDE = -I$(MSTOOLKIT) -I$(MSTOOLKIT)/mzParser
PEP = CKronik2.o


bullseye : bullseye.cpp $(PEP) 
	$(CC) $(FLAGS) $(INCLUDE) $(PEP) bullseye.cpp -L$(MSTOOLKIT) -lmstoolkitlite -o bullseye

#PEP objects
CKronik2.o : CKronik2.cpp
	$(CC) $(FLAGS) CKronik2.cpp -c

clean:
	rm -f *.o bullseye
