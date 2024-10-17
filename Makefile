CFLAGS=-O3 -std=c++11
#CFLAGS=-g

all: reducedlut.o
	g++ $(CFLAGS) -o reducedlut reducedlut.o

reducedlut.o: reducedlut.cpp reducedlut.h
	g++ $(CFLAGS) -c reducedlut.cpp

clean: 
	rm *.o reducedlut.exe
