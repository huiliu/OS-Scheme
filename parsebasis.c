/*
 * File: parsebasis.c
 * --------------------------------
 * This Program parse the basis function information which come from ESML.
 * Author:  liuhui
 * Date: Wed Sep  5 15:52:20 CST 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "parsebasis.h"
#include "print.h"
#include "common.h"

#define DEBUG_OUTPUT_INPUTFILE      1
#define DEBUG_INPUT                 0
#define DEBUG_OUTPUT_ALLBASIS       0
#define DEBUG_OUTPUT_ALLGTO         0
#define DEBUG_ZETA_K_P              0

#define INITIAL_ATOM_COUNT          100
#define INITIAL_BASIS_COUNT         200
#define INDEPENDENT_GTO_COUNT       300
#define PI_1_25         4.182513398379599

/*
    FUNCTION:  GET_BASIS_FILE(a, b)

    This code get the name of the Basis Function Database
    acoording the basis which user set
 */
inline void Get_Basis_File(char *basisName, char *BasisFile)
{
    if (strcasecmp(basisName, "6-31g*") == 0)
        strcat(BasisFile, "631gd");
    else if (strcasecmp(basisName, "6-31g") == 0)
        strcat(BasisFile, "631g");
    else if (strcasecmp(basisName, "6-31g**") == 0)
        strcat(BasisFile, "631gdp");
    else if (strcasecmp(basisName, "sto-3g") == 0)
        strcat(BasisFile, "sto3g");
    else if (strcasecmp(basisName, "acct") == 0)
        strcat(BasisFile, "acct");
    else if (strcasecmp(basisName, "cct") == 0)
        strcat(BasisFile, "cct");
}

INPUT_INFO* parse_input(const char* file_name)
{
/*
  hf 6-31g*

  0 1
  N    7        0.0000000000   0.0000000000   0.0000000000
  N    7        0.0000000000   0.0000000000   1.0980000000
*/
    int i = 0, j;
    FILE *f;
    INPUT_INFO *inputFile;
    char BasisFile[20] = "EMSL/";     // save the name of basis file
    char method[5], basisName[10];
    ATOM_INFO *atomList;
    COORD *coord;
    int gCount = 0, basisCount;


    // initial allocate 100 basis function, 300 GTO;
    inputFile = calloc(sizeof(INPUT_INFO), 1);
    atomList = inputFile->atomList =
                                   malloc(sizeof(ATOM_INFO)*INITIAL_ATOM_COUNT);
    coord = inputFile->gXYZ = malloc(sizeof(COORD) * INITIAL_ATOM_COUNT);

    OPEN_FILE(f, file_name);

    // read the control parameter like Gaussian 09
    fscanf(f, "%s%s", method, basisName);
    // Get the electronic count and multidegree of the system
    fscanf(f, "%d%u", &inputFile->icharge, &inputFile->imult);

#if DEBUG_OUTPUT_INPUTFILE
    fprintf(stdout, "%s %s\n\n%d %u\n",
                    method, basisName, inputFile->icharge, inputFile->imult);
#endif

    // read the atom information include atom Number and coordination
    while (fscanf(f, "%s%d%lf%lf%lf", atomList[i].symbol, &atomList[i].n,
                                &coord[i].x, &coord[i].y, &coord[i].z) != EOF) {
        // tranverse ANGS to BOHR
        coord[i].x *= ANGS_2_BOHR;
        coord[i].y *= ANGS_2_BOHR;
        coord[i].z *= ANGS_2_BOHR;
#if DEBUG_OUTPUT_INPUTFILE
        fprintf(stdout, "%s %d %lf %lf %lf\n",
                        atomList[i].symbol,
                        atomList[i].n,
                        coord[i].x,
                        coord[i].y,
                        coord[i].z);
#endif
        // connect the coordination
        inputFile->atomList[i].cid = i;
        inputFile->atomCount++;
        i++;
    }
    fclose(f);

    // release redundant memory
    REALLOC(inputFile->atomList, inputFile->atomCount * sizeof(ATOM_INFO));

    // ----------------------------------------------------------------
    // Read Basis information from the basis function database
    // ----------------------------------------------------------------
    inputFile->basisSet = malloc(sizeof(BASIS)* INITIAL_BASIS_COUNT);
    inputFile->gp = malloc(sizeof(GTO_PARTIAL)* INDEPENDENT_GTO_COUNT);
    inputFile->gtoSet = malloc(sizeof(GTO)* INDEPENDENT_GTO_COUNT);

    Get_Basis_File(basisName, BasisFile);
    OPEN_FILE(f, BasisFile);
    readbasis(f, inputFile);
    fclose(f);

    /* ----------------------------------------------------------------
     * allocate memory for store P, K, gamma
     * ----------------------------------------------------------------
     */
    gCount = inputFile->gCount;
    basisCount = inputFile->basisCount;
#if DEBUG_INPUT
    fprintf(stdout, "%d %d %d\n", gCount, inputFile->gtoCount,  basisCount);

    for (i = 0; i < gCount; i++)
        fprintf(stdout, "%d %lf %lf\n", i, inputFile->gp[i].alpha,
                                                        inputFile->gp[i].coeff);
#endif

    inputFile->K = (double **)malloc(sizeof(double *) * gCount);
    inputFile->zeta = (double **)malloc(sizeof(double *) * gCount);
    inputFile->P = (COORD **)malloc(sizeof(COORD *) * gCount);
    for (i = 0; i < gCount; i ++) {
        MALLOC(inputFile->K[i], sizeof(double) * gCount);
        MALLOC(inputFile->zeta[i], sizeof(double) * gCount);
        MALLOC(inputFile->P[i], sizeof(COORD) * gCount);
    }

    for (i = 0; i < gCount; i++) {
        for (j = 0; j <= i; j++) {
            // compute \zeta = \alpha_1 + \alpha_2
            inputFile->zeta[i][j] = inputFile->zeta[j][i] = 
                        inputFile->gp[i].alpha + inputFile->gp[j].alpha;
            // compute K = f(\alpha_1, \alpha_2, AB)
            inputFile->K[i][j] = inputFile->K[j][i] = K_OS(&inputFile->gp[i], 
                                                           &inputFile->gp[j],
                                                           inputFile->gXYZ);
            // compute P
            Gaussian_product_center(&inputFile->gp[i], &inputFile->gp[j],
                                                       inputFile->gXYZ,
                                                       &inputFile->P[i][j]);
            inputFile->P[j][i] = inputFile->P[i][j];

#if DEBUG_ZETA_K_P
            fprintf(stdout, "i, j: %d %d\nzeta = %lf\tK = %lf\tP: %lf %lf %lf\n",
                                                        i, j,
                                                        inputFile->zeta[i][j],
                                                        inputFile->K[i][j],
                                                        inputFile->P[i][j].x,
                                                        inputFile->P[i][j].y,
                                                        inputFile->P[i][j].z);
#endif
        }
    }
    return inputFile;
}

/*
 *  Function:   K_OS
 *
 *  The program compute the value of formula (47) in J. Chem. Phys. 84(7),3963
 *
 *  %Latex formula
 *  \begin{equation}
 *      K(\alpha_1, alpha_2, \mathbf{A}, \mathbf{B})
 *      = \frac{\sqrt{2}\pi^{\fact{5}{4}}}{\alpha_1 + \alpha_2}
 *        exp\left(-\frac{\alpha_1\alpha_2}{\alpha_1 + \alpha_2}
 *        \mathbf{AB}^2 \right)
 *  \end{equation}
 */
inline double K_OS(const GTO_PARTIAL *gp1,const GTO_PARTIAL *gp2,const COORD *c)
{
    int i = gp1->cid;
    int j = gp2->cid;
    double AB  = DISTANCE(c[i], c[j]);
    double zeta = gp1->alpha + gp2->alpha;

    return M_SQRT2 * PI_1_25 / zeta * exp(-gp1->alpha * gp2->alpha / zeta * AB);
}

/*
 *  Function:   readbasis
 *
 *  The program scan the basis data file for searching the basis information.
 */
#define LINE_LEN            80
#define DELIMITER           "****"
#define DELIMITER_LEN       4
#define ORBITAL_TYPE_COUNT  7

void readbasis(FILE * f, INPUT_INFO *inputFile)
{
/* 
 *  The structure of basis data file (EMSL gaussian 94)
****
H     0 
S   3   1.00                           
      3.42525091             0.15432897
      0.62391373             0.53532814
      0.16885540             0.44463454
****
He     0                               
S   3   1.00
      6.36242139             0.15432897
      1.15892300             0.53532814
      0.31364979             0.44463454
****       
Li     0       
*/
    int i, unkown;
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
                if (strcasecmp(symbol, inputFile->atomList[i].symbol) == 0)
                // the location of aim Atom was located.
                // Now start read basis information
                    GetBasis(f, inputFile, inputFile->atomList[i].cid);
            }
        }
                        rewind(f);
    }
#if DEBUG_OUTPUT_ALLBASIS
    PrintAllBasis(inputFile);
#endif
#if DEBUG_OUTPUT_ALLGTO
    gtoOutput(inputFile->gtoSet, inputFile->gXYZ, inputFile->gp,
                                                  inputFile->gtoCount);
#endif
}

/*
 *  Function:   GetBasis
 *
 *  The program read basis function parameter from the basis data file
 *  which come from EMSL.
 *
 */
int GetBasis(FILE *f, INPUT_INFO *inputFile, int cid)
{
/*
 *      S   3   1.00                           
 *      3.42525091             0.15432897
 *      0.62391373             0.53532814
 *      0.16885540             0.44463454
 */
    char symbol[5];
    char orbitalType[ORBITAL_TYPE_COUNT][3] =
                                        {"S", "P", "SP", "D", "F", "G", "H"};
    int state = 0, tmp_i = 0, i;
    double tmp_coeff = 0;
    GTO_PARTIAL *gp = inputFile->gp;
    GTO *g = inputFile->gtoSet;
    BASIS *b = inputFile->basisSet;
    int gCount = inputFile->gCount;
    int gtoCount = inputFile->gtoCount;
    int basisCount = inputFile->basisCount;

    while (1) {
        switch (state) {
            case 0: // initial state    S   3   1.00
                fscanf(f, "%s", symbol);
                if (strcmp(symbol, DELIMITER) == 0) {
                // the current basis block have done.
                    inputFile->gCount = gCount;
                    inputFile->gtoCount = gtoCount;
                    inputFile->basisCount = basisCount;
                    return 0;
                }

                fscanf(f, "%d%lf", &b[basisCount].gaussCount,
                                                        &b[basisCount].scale);
                // Get the orbital (basis function) type and goto next state
                for (i = 0; i < ORBITAL_TYPE_COUNT; i++)
                    if (strcmp(symbol, orbitalType[i]) == 0)
                        break;
                state = i + 1;
                break;
            case 1:
            // Read S orbital basis function 
                b[basisCount].gaussian = &g[gtoCount];
                b[basisCount].cid = cid;
                for (i = 0; i < b[basisCount].gaussCount; i++) {
                    fscanf(f, "%lf%lf", &gp[gCount].alpha, &gp[gCount].coeff);
                    gp[gCount].cid = cid;
                    g[gtoCount].gid = gCount;
                    g[gtoCount].coeff = Normalize(&g[gtoCount], &gp[gCount]);
                    gtoCount++;
                    gCount++;
                }
                // Use int type data indicate the orbital Type
                // strcpy(b[basisCount].type, "S");
                b[basisCount].Type = 0;
                basisCount++;
                state = 0;
                break;
            case 2:
            // Read P orbital basis function parameter
                tmp_i = b[basisCount].gaussCount;
                b[basisCount].gaussian = &g[gtoCount];
                b[basisCount + 1].gaussian = &g[gtoCount + tmp_i];
                b[basisCount + 2].gaussian = &g[gtoCount + tmp_i * 2];

                b[basisCount].cid =
                b[basisCount + 1].cid =
                b[basisCount + 2].cid = cid;

                b[basisCount + 1].gaussCount = 
                b[basisCount + 2].gaussCount = tmp_i;

                for (i = 0; i < tmp_i; i++) {
                    fscanf(f, "%lf%lf", &gp[gCount].alpha, &gp[gCount].coeff);
                    gp[gCount].cid = cid;

                    g[gtoCount].gid =
                    g[gtoCount + tmp_i].gid =
                    g[gtoCount + tmp_i * 2].gid = gCount;

                    g[gtoCount].l =
                    g[gtoCount + tmp_i].m = 
                    g[gtoCount + tmp_i * 2].n = 1;

                    g[gtoCount].coeff =
                    g[gtoCount + tmp_i].coeff = 
                    g[gtoCount + tmp_i * 2].coeff = Normalize(&g[gtoCount],
                                                              &gp[gCount]);

                    gtoCount++;
                    gCount++;
                }

                b[basisCount].Type = 1;
                b[basisCount + 1].Type = 2;
                b[basisCount + 2].Type = 3;

                gtoCount += tmp_i * 2;
                //gCount += tmp_i;
                basisCount += 3;

                state = 0;
                break;
            case 3:
            // Read SP orbital basis function parameter
                tmp_i = b[basisCount].gaussCount;

                b[basisCount].gaussian = &g[gtoCount];
                b[basisCount + 1].gaussian = &g[gtoCount + tmp_i];
                b[basisCount + 2].gaussian = &g[gtoCount + tmp_i * 2];
                b[basisCount + 3].gaussian = &g[gtoCount + tmp_i * 3];

                b[basisCount].cid =
                b[basisCount + 1].cid =
                b[basisCount + 2].cid =
                b[basisCount + 3].cid = cid;

                b[basisCount + 1].gaussCount = 
                b[basisCount + 2].gaussCount = 
                b[basisCount + 3].gaussCount = tmp_i;

                for (i = 0; i < tmp_i; i++) {
                    fscanf(f, "%lf%lf%lf", &gp[gCount].alpha, &gp[gCount].coeff,
                                                    &tmp_coeff);
                    // S
                    gp[gCount].cid = cid;
                    g[gtoCount].gid = gCount;
                    g[gtoCount].coeff = Normalize(&g[gtoCount], &gp[gCount]);

                    // P
                    gp[gCount + tmp_i].alpha = gp[gCount].alpha;
                    gp[gCount + tmp_i].coeff = tmp_coeff;
                    gp[gCount + tmp_i].cid = cid;

                    g[gtoCount + tmp_i].gid =
                    g[gtoCount + tmp_i * 2].gid =
                    g[gtoCount + tmp_i * 3].gid = gCount + tmp_i;


                    g[gtoCount + tmp_i].l = 
                    g[gtoCount + tmp_i * 2].m = 
                    g[gtoCount + tmp_i * 3].n = 1;

                    g[gtoCount + tmp_i].coeff =
                    g[gtoCount + tmp_i * 2].coeff = 
                    g[gtoCount + tmp_i * 3].coeff =
                           Normalize(&g[gtoCount + tmp_i], &gp[gCount + tmp_i]);

                    gtoCount++;
                    gCount++;
                }

                b[basisCount].Type = 0;
                b[basisCount + 1].Type = 1;
                b[basisCount + 2].Type = 2;
                b[basisCount + 3].Type = 3;

                gtoCount += tmp_i * 3;
                gCount += tmp_i;
                basisCount += 4;
                state = 0;
                break;
            case 4:
            // Read D orbital basis function parameter
                tmp_i = b[basisCount].gaussCount;

                b[basisCount].gaussian = &g[gtoCount];
                b[basisCount + 1].gaussian = &g[gtoCount + tmp_i];
                b[basisCount + 2].gaussian = &g[gtoCount + tmp_i * 2];
                b[basisCount + 3].gaussian = &g[gtoCount + tmp_i * 3];
                b[basisCount + 4].gaussian = &g[gtoCount + tmp_i * 4];
                b[basisCount + 5].gaussian = &g[gtoCount + tmp_i * 5];

                b[basisCount].cid =
                b[basisCount + 1].cid =
                b[basisCount + 2].cid =
                b[basisCount + 3].cid =
                b[basisCount + 4].cid =
                b[basisCount + 5].cid = cid;

                b[basisCount + 1].gaussCount =
                b[basisCount + 2].gaussCount =
                b[basisCount + 3].gaussCount =
                b[basisCount + 4].gaussCount =
                b[basisCount + 5].gaussCount = tmp_i;

                for (i = 0; i < tmp_i; i++) {
                    fscanf(f, "%lf%lf", &gp[gCount].alpha, &gp[gCount].coeff);
                    /*
                    g[gtoCount].gtoID = 
                    g[gtoCount + tmp_i].gtoID = 
                    g[gtoCount + tmp_i * 2].gtoID = 
                    g[gtoCount + tmp_i * 3].gtoID = 
                    g[gtoCount + tmp_i * 4].gtoID = 
                    g[gtoCount + tmp_i * 5].gtoID = (gtoCount) + tmp_i;
                    */
                    gp[gCount].cid =cid;

                    g[gtoCount].gid =
                    g[gtoCount + tmp_i].gid =
                    g[gtoCount + tmp_i * 2].gid =
                    g[gtoCount + tmp_i * 3].gid =
                    g[gtoCount + tmp_i * 4].gid =
                    g[gtoCount + tmp_i * 5].gid = gCount;

                    // d_{x^2} d_{y^2} d_{z^2}
                    g[gtoCount].l =
                    g[gtoCount + tmp_i].m = 
                    g[gtoCount + tmp_i * 2].n = 2;
                    // d_{xy} d_{xz} d_{yz}
                    g[gtoCount + tmp_i * 3].l =
                    g[gtoCount + tmp_i * 3].m =
                    g[gtoCount + tmp_i * 4].l =
                    g[gtoCount + tmp_i * 4].n =
                    g[gtoCount + tmp_i * 5].m =
                    g[gtoCount + tmp_i * 5].n = 1;

                    g[gtoCount].coeff =
                    g[gtoCount + tmp_i].coeff =
                    g[gtoCount + tmp_i * 2].coeff =
                                           Normalize(&g[gtoCount], &gp[gCount]);
                    g[gtoCount + tmp_i * 3].coeff =
                    g[gtoCount + tmp_i * 4].coeff =
                    g[gtoCount + tmp_i * 5].coeff =
                               Normalize(&g[gtoCount + tmp_i * 3], &gp[gCount]);

                    gtoCount++;
                    gCount++;
                }

                b[basisCount].Type = 4;
                b[basisCount + 1].Type = 5;
                b[basisCount + 2].Type = 6;
                b[basisCount + 3].Type = 7;
                b[basisCount + 4].Type = 8;
                b[basisCount + 5].Type = 9;

                gtoCount += tmp_i * 5;
                //gCount += tmp_i;
                basisCount += 6;
                state = 0;
                break;
            case 5:
            // Read F orbital basis function parameter
            // Use int type data indicate the orbital Type
                tmp_i = b[basisCount].gaussCount;

                b[basisCount].gaussian = &g[gtoCount];
                b[basisCount + 1].gaussian = &g[gtoCount + tmp_i];
                b[basisCount + 2].gaussian = &g[gtoCount + tmp_i * 2];
                b[basisCount + 3].gaussian = &g[gtoCount + tmp_i * 3];
                b[basisCount + 4].gaussian = &g[gtoCount + tmp_i * 4];
                b[basisCount + 5].gaussian = &g[gtoCount + tmp_i * 5];
                b[basisCount + 6].gaussian = &g[gtoCount + tmp_i * 6];
                b[basisCount + 7].gaussian = &g[gtoCount + tmp_i * 7];
                b[basisCount + 8].gaussian = &g[gtoCount + tmp_i * 8];
                b[basisCount + 9].gaussian = &g[gtoCount + tmp_i * 9];

                b[basisCount].cid =
                b[basisCount + 1].cid =
                b[basisCount + 2].cid =
                b[basisCount + 3].cid =
                b[basisCount + 4].cid =
                b[basisCount + 5].cid =
                b[basisCount + 6].cid =
                b[basisCount + 7].cid =
                b[basisCount + 8].cid =
                b[basisCount + 9].cid = cid;

                b[basisCount + 1].gaussCount =
                b[basisCount + 2].gaussCount =
                b[basisCount + 3].gaussCount =
                b[basisCount + 4].gaussCount =
                b[basisCount + 5].gaussCount =
                b[basisCount + 6].gaussCount =
                b[basisCount + 7].gaussCount =
                b[basisCount + 8].gaussCount =
                b[basisCount + 9].gaussCount = tmp_i;

                for (i = 0; i < tmp_i; i++) {
                    fscanf(f, "%lf%lf", &gp[gCount].alpha, &gp[gCount].coeff);

                    gp[gCount].cid =cid;

                    g[gtoCount].gid =
                    g[gtoCount + tmp_i].gid =
                    g[gtoCount + tmp_i * 2].gid =
                    g[gtoCount + tmp_i * 3].gid =
                    g[gtoCount + tmp_i * 4].gid =
                    g[gtoCount + tmp_i * 5].gid =
                    g[gtoCount + tmp_i * 6].gid =
                    g[gtoCount + tmp_i * 7].gid =
                    g[gtoCount + tmp_i * 8].gid =
                    g[gtoCount + tmp_i * 9].gid = gCount;

                    // f_{x^3} f_{y^3} f_{z^3}
                    g[gtoCount].l =
                    g[gtoCount + tmp_i].m = 
                    g[gtoCount + tmp_i * 2].n = 3;
                    // f_{x^2y}
                    g[gtoCount + tmp_i * 3].l = 2;
                    g[gtoCount + tmp_i * 3].m = 1;
                    //  f_{x^2z}
                    g[gtoCount + tmp_i * 4].l = 2;
                    g[gtoCount + tmp_i * 4].n = 1;
                    // f_{xy^2}
                    g[gtoCount + tmp_i * 5].l = 1;
                    g[gtoCount + tmp_i * 5].m = 2;
                    //  f_{y^2z}
                    g[gtoCount + tmp_i * 6].m = 2;
                    g[gtoCount + tmp_i * 6].n = 1;
                    // f_{xz^2}
                    g[gtoCount + tmp_i * 7].l = 1;
                    g[gtoCount + tmp_i * 7].n = 2;
                    //  f_{yz^2}
                    g[gtoCount + tmp_i * 8].m = 1;
                    g[gtoCount + tmp_i * 8].n = 2;
                    // f_{xyz}
                    g[gtoCount + tmp_i * 9].l = 1;
                    g[gtoCount + tmp_i * 9].m = 1;
                    g[gtoCount + tmp_i * 9].n = 1;

                    g[gtoCount].coeff =
                    g[gtoCount + tmp_i].coeff =
                    g[gtoCount + tmp_i * 2].coeff =
                                           Normalize(&g[gtoCount], &gp[gCount]);
                    g[gtoCount + tmp_i * 3].coeff =
                    g[gtoCount + tmp_i * 4].coeff =
                    g[gtoCount + tmp_i * 5].coeff =
                    g[gtoCount + tmp_i * 6].coeff =
                    g[gtoCount + tmp_i * 7].coeff =
                    g[gtoCount + tmp_i * 8].coeff =
                               Normalize(&g[gtoCount + tmp_i * 4], &gp[gCount]);
                    g[gtoCount + tmp_i * 9].coeff =
                               Normalize(&g[gtoCount + tmp_i * 9], &gp[gCount]);

                    gtoCount++;
                    gCount++;
                }

                b[basisCount].Type     = 10;
                b[basisCount + 1].Type = 11;
                b[basisCount + 2].Type = 12;
                b[basisCount + 3].Type = 13;
                b[basisCount + 4].Type = 14;
                b[basisCount + 5].Type = 15;
                b[basisCount + 6].Type = 16;
                b[basisCount + 7].Type = 17;
                b[basisCount + 8].Type = 18;
                b[basisCount + 9].Type = 19;

                gtoCount += tmp_i * 9;
                //gCount += tmp_i;
                basisCount += 10;

                state = 0;
                break;
            case 6:
            // Read G orbital basis function parameter
                state = 0;
                break;
            default:
                fprintf(stderr, "Read basis information error!\n");
                exit(EXIT_FAILURE);
        }
    }
}

/*
 *  Function:   Normalize
 *
 *  The program compute the normalization of the GTO.
 *
 *  Parameter:
 *  
 *      const GTO *g        contain the angular momentum(l, m, n).
 *      const GTO_PARTIAL *gp       contain the exponent and the contract
 *                                  coefficient
 */
inline double Normalize(const GTO *g, const GTO_PARTIAL *gp)
{
    double alpha = gp->alpha;
    int l = g->l;
    int m = g->m;
    int n = g->n;

    return gp->coeff * pow(2*alpha/M_PI, 0.75) * sqrt(pow(4*alpha, l + m + n)/
        (factorial_2(2*l-1) * factorial_2(2*m-1) * factorial_2(2*n-1)));
}
