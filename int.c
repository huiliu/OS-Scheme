#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "print.h"
#include "common.h"
#include "parsebasis.h"
#include "eri_drive.h"

int main(int argc, char** argv)
{
    char *inp = NULL;
    clock_t start, finish;
    int n;

    if (argc < 2)
        inp = "new.inp";
    else
        inp = argv[1];

    INPUT_INFO *b = parse_input(inp);    
    n = b->basisCount;
    //printf("Basis Count: %d\n", n);

/*
    gsl_matrix* overlapMatrix = overlap_matrix(b);
    matrix_output(overlapMatrix, n, "OVERLAP INTEGRALS:");

    gsl_matrix* kmatrix = kinetic_matrix(b);
    matrix_output(kmatrix, n, "KINETIC INTEGRALS:");

    gsl_matrix* attractionmatrix = nuclear_attraction_matrix(b);
    matrix_output(attractionmatrix, n, "ATTRACTION INTEGRALS:");

    gsl_matrix* h = hamiltonian(b);
    matrix_output(h, n, "Hamiltonian Matirx:");
*/
    double ****int2e;
    start = clock();
    int2e = ERI_Matrix(b);
    finish = clock();
    fprintf(stderr, "ERI CPU Time: %f\n", (double)(finish - start) / CLOCKS_PER_SEC);
    int2e_output(int2e, n, "TWO ELECTRON INTEGRAL:");

    return 0;
}
