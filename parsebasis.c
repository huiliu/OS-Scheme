/*
 * File: parsebasis.c
 * --------------------------------
 * This Program parse the basis function information which come from ESML.
 * Author:  liuhui
 * Date: Wed Sep  5 15:52:20 CST 2012
 */

/*
  hf 6-31g*

  0 1
  N    7.0        0.0000000000   0.0000000000   0.0000000000
  N    7.0        0.0000000000   0.0000000000   1.0980000000
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "parsebasis.h"
#include "print.h"

#define OPEN_FILE(f,fName)  \
    if((f = fopen(fName, "r")) == NULL) { \
        fprintf(stdout, "Open file \"%s\" failed!\n", fName); \
        exit(EXIT_FAILURE); \
    }

/*
    FUNCTION:  GET_BASIS_FILE(a, b)
    This code get the file name acoording the basis which user set
 */
inline void Get_Basis_File(char *basisName, char *BasisFile)
{
    if (strcmp(basisName, "6-31g*") == 0)
        strcat(BasisFile, "631gd");
    else if (strcmp(basisName, "sto-3g") == 0)
        strcat(BasisFile, "sto3g");
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stdout, "USAGE:\n   ./pasrebasis <basisFile>\n");
        exit(EXIT_FAILURE);
    }

    parse_input(argv[1]);

    return 0;
}

#define INITIAL_ATOM_COUNT          100
#define INITIAL_BASIS_COUNT         100
#define INITIAL_GTO_COUNT           300

INPUT_INFO* parse_input(const char* file_name)
{
    int i = 0;
    FILE *f;
    INPUT_INFO *inputFile;
    char BasisFile[20] = "EMSL/";     // save the name of basis file
    char method[5], basisName[10];
    ATOM_INFO *atomList;
    COORD *coord;


    // initial allocate 100 basis function, 300 GTO;
    inputFile = calloc(sizeof(INPUT_INFO), 1);
    atomList = inputFile->atomList = malloc(sizeof(ATOM_INFO)*INITIAL_ATOM_COUNT);
    inputFile->basisSet = malloc(sizeof(BASIS)* INITIAL_BASIS_COUNT);
    inputFile->gtoSet = malloc(sizeof(BASIS)* INITIAL_GTO_COUNT);
    coord = inputFile->gXYZ = malloc(sizeof(COORD)* INITIAL_ATOM_COUNT);
    inputFile->P = (COORD **)malloc(sizeof(COORD *)* INITIAL_GTO_COUNT);

    OPEN_FILE(f, file_name);

    // read the control parameter like Gaussian 09
    fscanf(f, "%s%s", method, basisName);
    fscanf(f, "%d%u", &inputFile->icharge, &inputFile->imult);

    // read the atom information include atom Number and coordination
    while (fscanf(f, "%s%d%lf%lf%lf", atomList[i].symbol, &atomList[i].n,
                                &coord[i].x, &coord[i].y, &coord[i].z) != EOF) {
        fprintf(stdout, "%s %s %d %d %s %d %lf %lf %lf\n",
                        method, basisName,
                        inputFile->icharge, inputFile->imult,
                        atomList[i].symbol, atomList[i].n,
                        coord[i].x, coord[i].y, coord[i].z);
        inputFile->atomCount++;
        i++;
    }
    fclose(f);

    // ----------------------------------------------------------------
    // Read Basis information from the basis function database
    Get_Basis_File(basisName, BasisFile);
    OPEN_FILE(f, BasisFile);
    readbasis(f, inputFile);
    fclose(f);
    return NULL;
}
    
#define LINE_LEN            80
#define DELIMITER           "****"
#define DELIMITER_LEN       4
#define ORBITAL_TYPE_COUNT  6

//BASIS* readbasis(FILE * f, BASIS *b, int basisCount, GTO *g, int *gtoCount)
void readbasis(FILE * f, INPUT_INFO *inputFile)
{
    int i, unkown, gtoCount = 0, basisCount = 0;
    char string[81], symbol[3];

    for (i = 0; i < inputFile->atomCount; i++) {
        while(fgets(string, LINE_LEN, f)) {
            if (feof(f)) break;
            // ignore the information until find the first DELIMITER
            if (strncmp(string, DELIMITER, DELIMITER_LEN) == 0) {
                // read the elemental symbol
                //fgets(string, LINE_LEN, f);
                fscanf(f, "%s%d", symbol, &unkown);

                /* TODO:
                 *      Optimize 
                 */
                    if (strcasecmp(symbol, inputFile->atomList[i].symbol) == 0) {
                    // the location of aim Atom was located.
                    // Now start read basis information
                        GetBasis(f, inputFile->gtoSet, &gtoCount, 
                                    inputFile->basisSet, &basisCount);
                    }
            }
        }
                        rewind(f);
    }
    gtoOutput(inputFile->gtoSet, gtoCount);
}
 
int GetBasis(FILE *f, GTO *g, int *gStart, BASIS *b, int *bStart)
{
    char symbol[5];
    char orbitalType[ORBITAL_TYPE_COUNT][3] = {"S", "SP", "D", "F", "G", "H"};
    int state = 0, tmp_i = 0, i;
/*
 *      S   3   1.00                           
 *      3.42525091             0.15432897
 *      0.62391373             0.53532814
 *      0.16885540             0.44463454
 */
    while (1) {
        switch (state) {
            case 0: // initial state    S   3   1.00
                fscanf(f, "%s", symbol);
                if (strcmp(symbol, DELIMITER) == 0)
                // the current basis block have done.
                    return 0;

                strcpy(b[*bStart].type, symbol);
                fscanf(f, "%d%lf", &b[*bStart].gaussCount, &b[*bStart].scale);
                // Get the orbital (basis function) type and goto next state
                for (i = 0; i < ORBITAL_TYPE_COUNT; i++)
                    if (strcmp(symbol, orbitalType[i]) == 0)
                        break;
                state = i + 1;
                break;
            case 1:
            // Read S orbital basis function 
                b[*bStart].gaussian = &g[*gStart];
                for (i = 0; i < b[*bStart].gaussCount; i++) {
                    fscanf(f, "%lf%lf", &g[*gStart].alpha, &g[*gStart].coeff);
                    g[*gStart].gtoID = *gStart;
                    g[*gStart].coeff *= Normalize(&g[*gStart]);
                    (*gStart)++;
                }
                (*bStart)++;
                state = 0;
                break;
            case 2:
            // Read SP orbital basis function parameter
                tmp_i = b[*bStart].gaussCount;

                b[*bStart].gaussian = &g[*gStart];
                b[*bStart + 1].gaussian = &g[*gStart + tmp_i];
                b[*bStart + 2].gaussian = &g[*gStart + tmp_i * 2];
                b[*bStart + 3].gaussian = &g[*gStart + tmp_i * 3];

                for (i = 0; i < tmp_i; i++) {
                    fscanf(f, "%lf%lf%lf", &g[*gStart].alpha, &g[*gStart].coeff,
                                                    &g[*gStart+tmp_i].coeff);
                    // S
                    g[*gStart].gtoID = *gStart;
                    g[*gStart].coeff *= Normalize(&g[*gStart]);

                    // P
                    g[*gStart + tmp_i].gtoID = 
                    g[*gStart + tmp_i * 2].gtoID = 
                    g[*gStart + tmp_i * 3].gtoID = (*gStart) + tmp_i;

                    g[*gStart + tmp_i].alpha = 
                    g[*gStart + tmp_i * 2].alpha = 
                    g[*gStart + tmp_i * 3].alpha = g[*gStart].alpha;

                    g[*gStart + tmp_i].l = 
                    g[*gStart + tmp_i * 2].m = 
                    g[*gStart + tmp_i * 3].n = 1;

                    g[*gStart + tmp_i].coeff *= Normalize(&g[*gStart + tmp_i]);
                    g[*gStart + tmp_i * 2].coeff = 
                    g[*gStart + tmp_i * 3].coeff = g[*gStart + tmp_i].coeff;

                    (*gStart)++;
                }
                (*gStart) += tmp_i * 3;
                (*bStart) += 4;
                state = 0;
                break;
            case 3:
            // Read D orbital basis function parameter
                tmp_i = b[*bStart].gaussCount;

                b[*bStart].gaussian = &g[*gStart];
                b[*bStart + 1].gaussian = &g[*gStart + tmp_i];
                b[*bStart + 2].gaussian = &g[*gStart + tmp_i * 2];
                b[*bStart + 3].gaussian = &g[*gStart + tmp_i * 3];
                b[*bStart + 4].gaussian = &g[*gStart + tmp_i * 4];
                b[*bStart + 5].gaussian = &g[*gStart + tmp_i * 5];

                for (i = 0; i < tmp_i; i++) {
                    fscanf(f, "%lf%lf", &g[*gStart].alpha, &g[*gStart].coeff);

                    g[*gStart].gtoID = 
                    g[*gStart + tmp_i].gtoID = 
                    g[*gStart + tmp_i * 2].gtoID = 
                    g[*gStart + tmp_i * 3].gtoID = 
                    g[*gStart + tmp_i * 4].gtoID = 
                    g[*gStart + tmp_i * 5].gtoID = (*gStart) + tmp_i;

                    g[*gStart + tmp_i].alpha = 
                    g[*gStart + tmp_i * 2].alpha = 
                    g[*gStart + tmp_i * 3].alpha = 
                    g[*gStart + tmp_i * 4].alpha = 
                    g[*gStart + tmp_i * 5].alpha = g[*gStart].alpha;

                    // d_{x^2} d_{y^2} d_{z^2}
                    g[*gStart].l =
                    g[*gStart + tmp_i].m = 
                    g[*gStart + tmp_i * 2].n = 2;
                    // d_{xy} d_{xz} d_{yz}
                    g[*gStart + tmp_i * 3].l =
                    g[*gStart + tmp_i * 3].m =
                    g[*gStart + tmp_i * 4].l =
                    g[*gStart + tmp_i * 4].n =
                    g[*gStart + tmp_i * 5].m =
                    g[*gStart + tmp_i * 5].n = 1;

                    g[*gStart + tmp_i * 3].coeff = g[*gStart].coeff;
                    g[*gStart + tmp_i * 3].coeff *= 
                                                Normalize(&g[*gStart + tmp_i * 3]);
                    g[*gStart].coeff *= Normalize(&g[*gStart]);
                    g[*gStart + tmp_i].coeff =
                    g[*gStart + tmp_i * 2].coeff = g[*gStart].coeff;
                    g[*gStart + tmp_i * 4].coeff =
                    g[*gStart + tmp_i * 5].coeff = g[*gStart + tmp_i * 3].coeff;
                    (*gStart)++;
                }
                (*gStart) += tmp_i * 5;
                (*bStart) += 6;
                state = 0;
                break;
            case 4:
            // Read F orbital basis function parameter
                state = 0;
                break;
            case 5:
            // Read G orbital basis function parameter
                state = 0;
                break;
            default:
                fprintf(stderr, "Read basis information error!\n");
                exit(EXIT_FAILURE);
        }
    }
}

inline double Normalize(const GTO *g)
{
    double fact2[] = {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025, 10321920, 34459425, 185794560, 654729075, 3715891200, 13749310575, 81749606400, 316234143225, 1961990553600};
    double alpha = g->alpha;
    int l = g->l;
    int m = g->m;
    int n = g->n;

    if (l <= 0) l = 1;
    if (m <= 0) m = 1;
    if (n <= 0) n = 1;
    return pow(2 * alpha / M_PI, 0.75) * sqrt(pow(4*alpha, l + m + n) / \
        (fact2[2*l-1] * fact2[2*m-1] * fact2[2*n-1]));
}
