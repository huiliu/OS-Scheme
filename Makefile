cc = gcc
OPT = -g -Wall -lm -o

all: parsebasis print.o common.o

parsebasis: parsebasis.c parsebasis.h print.o common.o
	$(cc) parsebasis.c print.o common.o $(OPT) parsebasis

print.o: print.c print.h parsebasis.h
	$(cc) print.c -c $(OPT) print.o

common.o: common.c common.h
	$(cc) common.c -c $(OPT) common.o

clean:
	rm -rf parsebasis *.o
