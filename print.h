/*
 * File: print.c
 * ----------------------------------------
 * This file contain output function.
 *
 * Author: liuhui
 * Date: Sat Sep  8 00:12:46 CST 2012
 */

#ifndef ERI__PRINT__
#define ERI__PRINT__

#include "common.h"
void gtoOutput(const GTO *g, const COORD *c, const GTO_PARTIAL *gp, int n);
void int2e_output(double**** e, int n, char* msg);
void PrintAllBasis(const INPUT_INFO *inp);
void PrintBasis(const BASIS *b, const GTO_PARTIAL *gp, const COORD *coord);

#endif
