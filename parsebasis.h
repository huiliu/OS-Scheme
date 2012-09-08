/*
 * File: parsebasis.h
 * --------------------------------
 * This Program parse the basis function information which come from ESML.
 * Author:  liuhui
 * Date: Wed Sep  5 15:52:20 CST 2012
 */

#include "common.h"


#ifndef __INTEGRAL__BASIS__
#define __INTEGRAL__BASIS__

#define OPEN_FILE(f,fName)  \
    if((f = fopen(fName, "r")) == NULL) { \
        fprintf(stdout, "Open file \"%s\" failed!\n", fName); \
        exit(EXIT_FAILURE); \
    }

INPUT_INFO* parse_input(const char* file_name);
void readbasis(FILE * f, INPUT_INFO *inputFile);
int GetBasis(FILE *f, INPUT_INFO *inputFile, int cid);
inline double Normalize(const GTO *g, const GTO_PARTIAL *gp);
inline void Get_Basis_File(char *basisName, char *BasisFile);
inline double K_OS(const GTO_PARTIAL *gp1, const GTO_PARTIAL *gp2, const COORD *c);


#endif
