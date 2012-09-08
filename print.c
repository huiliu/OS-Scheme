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
#include "parsebasis.h"

void gtoOutput(const GTO *g, int n)
{
    int i;

    for (i = 0; i < n; i++)
    fprintf(stdout, "gtoID: %d\n %d %d %d %lf %lf\n\n", g[i].gtoID,
                g[i].l, g[i].m, g[i].n, g[i].alpha, g[i].coeff);
}
