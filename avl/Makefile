cc = gcc

OPT = -g -Wall -o

all: test example2 avl.o redblack.o

test: avl.o test.c
	$(cc) test.c avl.o $(OPT) test

example2: example2.c redblack.o
	$(cc) example2.c redblack.o $(OPT) example2

avl.o: avl.c avl.h
	$(cc) avl.c -c $(OPT) avl.o

redblack.o: redblack.c redblack.h
	$(cc) redblack.c -c $(OPT) redblack.o

clean:
	rm -rf avl.o test redblack.o example2
