cc = gcc
OPT = -g -Wall -lm -o

all: parsebasis print.o

parsebasis: parsebasis.c parsebasis.h print.o
	$(cc) parsebasis.c print.o $(OPT) parsebasis

print.o: print.c print.h parsebasis.h
	$(cc) print.c -c $(OPT) print.o

clean:
	rm -rf parsebasis *.o
