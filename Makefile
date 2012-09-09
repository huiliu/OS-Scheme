cc = gcc
OPT = -g -Wall -lm -o

all:  int

parsebasis.o: parsebasis.c parsebasis.h
	$(cc) parsebasis.c print.o common.o -c $(OPT) parsebasis.o

print.o: print.c print.h parsebasis.h
	$(cc) print.c -c $(OPT) print.o

common.o: common.c common.h
	$(cc) common.c -c $(OPT) common.o

int: int.c print.o common.o parsebasis.o
	$(cc) int.c print.o common.o parsebasis.o $(OPT) int

clean:
	rm -rf int *.o
