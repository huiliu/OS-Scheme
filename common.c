#include "common.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/*
 *  放置一些常用函数
 */

/*
 * Constant Value
 * ------------------------------------------
 */

double fact2[] = {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025, 10321920, 34459425, 185794560, 654729075, 3715891200, 13749310575, 81749606400, 316234143225, 1961990553600};

#define FACTORIAL_2(n)  ((n) <= 0 ? 1 : fact2[n])

#define M_SQRT_PI_2             1.253314137315500121
#define F_INC_GAMMA_CYCLE       100
#define F_INC_GAMMA_delta       1.0E-10
/*
 *  Function:   F_inc_gamma
 *  Compute the value of Boys function.
 *  This program's efficency is very important.
 *
 *  %Latex Math Formula
 *  \begin{equation}
 *      F_m(w) = \int_0^1 exp(-wt^2) t^{2m} dt
 *  \end{equation}
 */
inline double F_inc_gamma(int m ,double w)
{
    double result = 0;
    double tmp = 0;
    int i;
    
    if (w < 17) {
        result = tmp = 1.0 / FACTORIAL_2(2*m + 1);
        for (i = 1; i < F_INC_GAMMA_CYCLE; i++) {
            tmp *= ((2*w) / (2*m + 2*i + 1));
            if ((tmp - F_INC_GAMMA_delta) < 0)
                break;
            result += tmp;
        }
        return result * FACTORIAL_2(2 * m -1) * exp(-w);;
    }else
        result = FACTORIAL_2(2*m -1) / pow(2*w, m + 0.5) * M_SQRT_PI_2;
    return result;
}

inline int factorial(int n)
{
    int i, result = 1;
    if (n <= 1) return 1;
    for (i = 2; i <= n; i++)
        result *= i;
    return result;
}


inline int factorial_2(int n)
{
    if (n <= 0) return 1;
    return fact2[n];
}


// check the symtery of two-electron integral
#define GetIndex(i, j, k, l, N) (N * N * N * i + N * N * j + N * k + l)
#ifndef MIN
#define MIN(a, b)   ((a) < (b) ? (a) : (b))
#endif

/*
 *  Function:   ChkERISym
 *  This program check the symmetry of ERI
 *
 *  e           a four degree pointer store the ERI.
 *  i           the index of basis function
 *  j           the index of basis function
 *  k           the index of basis function
 *  l           the index of basis function
 *  N           dimension of ERI
 *  is_dup      return value indicate the ERI whether have done.
 *
 */
inline int ChkERISym(double ****e, int i, int j, int k, int l, int N, int *is_dup)
{
// thanks David pulq for his idea.
// https://plus.google.com/106075773891428215861/posts
// https://gist.github.com/3427265
    long pos;
    long min;

    if (e[i][j][k][l] != 0){
        *is_dup = 1;
        return 0;
    }

    pos = min = GetIndex(i, j, k, l, N);
    
    min = MIN(min, GetIndex(i, j, l, k, N));
    min = MIN(min, GetIndex(j, i, l, k, N));
    min = MIN(min, GetIndex(j, i, k, l, N));
    min = MIN(min, GetIndex(k, l, i, j, N));
    min = MIN(min, GetIndex(k, l, j, i, N));
    min = MIN(min, GetIndex(l, k, i, j, N));
    min = MIN(min, GetIndex(l, k, j, i, N));
    
    *is_dup = pos != min;
    
    return 0;
}

/*
 * Function:    Gaussian_product_center
 * This function compute the the center of the production of two gaussian function
 *
 *  a           the exponent of the first gaussian function.
 *  A           the coordination of the first gaussian function.
 *  b           the exponent of the first gaussian function.
 *  B           the coordination of the first gaussian function.
 *  p           the center of production.
 */
inline void Gaussian_product_center(const GTO_PARTIAL *gpa, const GTO_PARTIAL *gpb,  const COORD *c, COORD *p)
{
// Gaussian函数乘积定理计算双中心
    double alpha1 = gpa->alpha;
    double alpha2 = gpb->alpha;
    double zeta = alpha1 + alpha2;
    int i = gpa->cid;
    int j = gpb->cid;
    
    p->x = (alpha1 * c[i].x + alpha2 * c[j].x) / zeta;
    p->y = (alpha1 * c[i].y + alpha2 * c[j].y) / zeta;
    p->z = (alpha1 * c[i].z + alpha2 * c[j].z) / zeta;
}

inline int GetAngularMomentum(int i)
{
    if (i == 0) //  S
        return 0;
    else if (i <= 3 && i >= 1)  // P
        return 1;
    else if (i <= 9 && i >= 4)  // D
        return 2;
    else if (i <= 19 && i >=10) // F
        return 3;
    else {
        fprintf(stderr, "\"%d\" The type of orbital doesn't exsit.\n", i);
        exit(EXIT_FAILURE);
    }
}
