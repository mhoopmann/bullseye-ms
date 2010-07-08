#Set these variables if needed
C = gcc
CC = g++
FLAGS = -O3 -static -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILEOFFSET_BITS=64
MSTOOLKIT = ../MSToolkit
INCLUDE = -I$(MSTOOLKIT) -I$(MSTOOLKIT)/zLib -I$(MSTOOLKIT)/RAMP
PEP = cPersistPep.o


bullseye : bullseye.cpp $(PEP) 
	$(CC) $(FLAGS) $(INCLUDE) $(PEP) bullseye.cpp -L$(MSTOOLKIT) -lmstoolkitlite -o bullseye

#PEP objects
cPersistPep.o : cPersistPep.cpp
	$(CC) $(FLAGS) cPersistPep.cpp -c

clean:
	rm -f *.o bullseye
