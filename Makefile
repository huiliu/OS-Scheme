cc = icc
OPT = -g -pg -Wall -lm -o
#OPT = -O3 -fast -lm -o
#OPT = -O3 -pg -lm -o

all:  int

parsebasis.o: parsebasis.c parsebasis.h
	$(cc) parsebasis.c -c $(OPT) parsebasis.o

print.o: print.c print.h parsebasis.h
	$(cc) print.c -c $(OPT) print.o

common.o: common.c common.h
	$(cc) common.c -c $(OPT) common.o

eri_os.o: eri_os.c eri_os.h common.h print.h
	$(cc) eri_os.c -c $(OPT) eri_os.o

eri_hgp.o: eri_hgp.c eri_hgp.h common.h print.h
	$(cc) eri_hgp.c -c $(OPT) eri_hgp.o

eri_drive.o: eri_drive.c eri_drive.h common.h
	$(cc) eri_drive.c -c $(OPT) eri_drive.o

int: int.c print.o common.o parsebasis.o eri_drive.o eri_os.o eri_hgp.o
	$(cc) int.c print.o common.o parsebasis.o eri_drive.o eri_os.o eri_hgp.o $(OPT) int

clean:
	rm -rf int *.o gmon.out
