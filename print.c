/*
 * File: print.c
 * ----------------------------------------
 * This file contain output function.
 *
 * Author: liuhui
 * Date: Sat Sep  8 00:12:46 CST 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

void gtoOutput(const GTO *g, const COORD *c, const GTO_PARTIAL *gp, int n)
{
    int i;

    for (i = 0; i < n; i++)
    fprintf(stdout, "%d %d %d\n \
    gid: %d\n %15.9E %15.9E\n \
    cid: %d\n %15.9E %15.9E %15.9E\n\n", 
                g[i].l, g[i].m, g[i].n,
                g[i].gid, gp[g[i].gid].alpha, gp[g[i].gid].coeff, 
                gp[g[i].gid].cid,
                c[gp[g[i].gid].cid].x,
                c[gp[g[i].gid].cid].y,
                c[gp[g[i].gid].cid].z); 
}
void int2e_output(double**** e, int n, char* msg)
{
    int i, j, k, l;

    printf("%s\n", msg);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
                    if (fabs(e[i][j][k][l]) > 1.0E-8)
                        printf("%3d%3d%3d%3d%15.9lf\n", i, j, k, l, e[i][j][k][l]);
                }
            }
        }
    }
}

void PrintAllBasis(const INPUT_INFO *inp)
{
    int i, j, gaussCount, basisCount = inp->basisCount;
    BASIS *b = inp->basisSet;
    GTO_PARTIAL *gp = inp->gp;
    GTO *g;
    
    for (i = 0; i < basisCount; i++) {
        fprintf(stdout, "No. %d\nOrbital Type: %d\nThe count of Gaussian: %d\n",
                                                i, b[i].Type, b[i].gaussCount);
        g = b[i].gaussian;
        fprintf(stdout, "The coordination is:\n %d %16.9E %16.9E %16.9E\n", b[i].cid,
                    inp->gXYZ[b[i].cid].x, inp->gXYZ[b[i].cid].y,
                    inp->gXYZ[b[i].cid].z); 
        gaussCount = b[i].gaussCount;
        for (j = 0; j < gaussCount; j++) {
            fprintf(stdout, "%d %d %d %16.9E %16.9E %16.9E\n", g[j].l, g[j].m, g[j].n,
                        gp[g[j].gid].alpha, g[j].coeff, gp[g[j].gid].coeff);
        }

    }
}

void PrintBasis(const BASIS *b, const GTO_PARTIAL *gp, const COORD *coord)
{
    int i;
    int gaussCount = b->gaussCount;

    fprintf(stdout, "基组坐标:\n%15.06lE%15.06lE%15.06lE\n基函数\n",
                            coord[b->cid].x, coord[b->cid].y, coord[b->cid].z);
    for (i = 0; i < gaussCount; i++) {
        fprintf(stdout, "%d %d %d %12.8lf %12.8lf %12.8lf\n",
                    b->gaussian[i].l, b->gaussian[i].m, b->gaussian[i].n,
                    gp[b->gaussian[i].gid].alpha, b->gaussian[i].coeff,
                    gp[b->gaussian[i].gid].coeff);
    }
}
