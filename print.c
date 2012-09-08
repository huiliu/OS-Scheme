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
#include "common.h"

void gtoOutput(const GTO *g, const COORD *c, const GTO_PARTIAL *gp, int n)
{
    int i;

    for (i = 0; i < n; i++)
    fprintf(stdout, "%d %d %d\n \
    gid: %d\n %lf %lf\n \
    cid: %d\n %lf %lf %lf\n\n", 
                g[i].l, g[i].m, g[i].n,
                g[i].gid, gp[g[i].gid].alpha, gp[g[i].gid].coeff, 
                gp[g[i].gid].cid,
                c[gp[g[i].gid].cid].x,
                c[gp[g[i].gid].cid].y,
                c[gp[g[i].gid].cid].z); 
}
