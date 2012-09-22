#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eri_hgp.h"


double ERI_HGP_DRIVE(const INPUT_INFO *inp, int ii, int jj, int kk, int ll, int debug)
{
/*
 * select all basis function in common shell. for example px, py, pz;
 *
 */
    double *XSXS = NULL;
    double result;
    COORD AB, CD;
    BASIS b1, b2, b3, b4;

    b1 = inp->basisSet[ii];
    b2 = inp->basisSet[jj];
    b3 = inp->basisSet[kk];
    b4 = inp->basisSet[ll];

    PP(inp->gXYZ[b1.cid], inp->gXYZ[b2.cid], AB);
    PP(inp->gXYZ[b3.cid], inp->gXYZ[b4.cid], CD);

if (debug == 9) {
    fprintf(stdout, "AB:\t%12.6E%12.6E%12.6E\n", AB.x, AB.y, AB.z);
    fprintf(stdout, "CD:\t%12.6E%12.6E%12.6E\n", CD.x, CD.y, CD.z);
}

    result = HGPBasisHRR(inp, ii, jj, kk, ll, 
                             b1.gaussian[0].l, b1.gaussian[0].m, b1.gaussian[0].n,
                             b2.gaussian[0].l, b2.gaussian[0].m, b2.gaussian[0].n,
                             b3.gaussian[0].l, b3.gaussian[0].m, b3.gaussian[0].n,
                             b4.gaussian[0].l, b4.gaussian[0].m, b4.gaussian[0].n,
                             &AB, &CD, XSXS, debug);
                        
    return result;
    //free(XSXS);
}

#define COORDINATION_THRESHOLD          1.0E-12
double HGPBasisHRR(const INPUT_INFO *inp,
                    int ii, int jj, int kk, int ll, 
                    int l1, int m1, int n1,
                    int l2, int m2, int n2,
                    int l3, int m3, int n3,
                    int l4, int m4, int n4,
                    COORD *AB, COORD *CD, double *XSXS, int debug)
{
/* the input basis may be:
 * (px, py|dx^2, dxy), ......
 * transpose it to the form of (e0|f0) by HRR(J. Chem. Phys. 89(9), 5777, 1988)
 * formula (18)
 */

// contract basis to form (e0|f0)
    double item = 0;

// (ab|c(d+1)) = (ab|(c+1)d) + CD(ab|cd)
    if (l4 > 0) {
        item = HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2,n2,l3+1,m3,n3,l4-1,m4,n4,AB, CD, XSXS, debug);
        if (fabs(CD->x) > COORDINATION_THRESHOLD)
            item += CD->x * HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2,n2,l3,m3,n3,l4-1,m4,n4,AB, CD, XSXS, debug);
        return item;   
    }
    if (m4 > 0) {
        item = HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2,n2,l3,m3+1,n3,l4,m4-1,n4,AB, CD, XSXS, debug);
        if (fabs(CD->y) > COORDINATION_THRESHOLD)
            item += CD->y * HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2,n2,l3,m3,n3,l4,m4-1,n4,AB, CD, XSXS, debug);
        return item;
    }
    if (n4 > 0) {
        item = HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2,n2,l3,m3,n3+1,l4,m4,n4-1,AB, CD, XSXS, debug);
        if (fabs(CD->z) > COORDINATION_THRESHOLD)
            item += CD->z * HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2,n2,l3,m3,n3,l4,m4,n4-1,AB, CD, XSXS, debug);
        return item;
    }

 // (a(b+1)|cd) = ((a+1)b|cd) + AB(ab|cd)
    if (l2 > 0) {
        item = HGPBasisHRR(inp, ii, jj, kk, ll, l1+1,m1,n1,l2-1,m2,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        if (fabs(AB->x) > COORDINATION_THRESHOLD)
            item += AB->x * HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2-1,m2,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        return item;
    }
    if (m2 > 0) {
        item = HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1+1,n1,l2,m2-1,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        if (fabs(AB->y) > COORDINATION_THRESHOLD)
        item += AB->y * HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2-1,n2,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        return item;   
    }
    if (n2 > 0) {
        item = HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1+1,l2,m2,n2-1,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        if (fabs(AB->z) > COORDINATION_THRESHOLD)
            item += AB->z * HGPBasisHRR(inp, ii, jj, kk, ll, l1,m1,n1,l2,m2,n2-1,l3,m3,n3,l4,m4,n4, AB, CD, XSXS, debug);
        return item;
    }

    // the basis integral has been converted to the form of (e0|f0), 
    // and continue to carry out transpose (e0|f0) to [e0|f0]
    // FORM (e0|f0) FROM ∑∑∑∑[e0|f0]
    return HGPBasis(inp, ii, jj, kk, ll, l1, m1, n1, l3, m3, n3, XSXS, debug);
}

double HGPBasis(const INPUT_INFO *inp,
                int ii, int jj, int kk, int ll,
                int l1, int m1, int n1,
                int l3, int m3, int n3,
                double *XSXS, int debug)
{
    // [e0|f0] ---> (e0|f0)
    int i, j, k, l, gaussCount_1, gaussCount_2, gaussCount_3, gaussCount_4;
    int gid1, gid2, gid3, gid4;
    int cid1, cid2, cid3, cid4;
    int L, m;
    double zeta, eta, rho;
    BASIS *b1, *b2, *b3, *b4;
    GTO *g1, *g2, *g3, *g4;
    COORD PA, PB, QC, QD, WQ, WP, PQ; 
    double KAB, KCD, pre;
    double *F;
    double T;
    double result = 0;

    b1 = &inp->basisSet[ii];
    b2 = &inp->basisSet[jj];
    b3 = &inp->basisSet[kk];
    b4 = &inp->basisSet[ll];

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;
    gaussCount_3 = b3->gaussCount;
    gaussCount_4 = b4->gaussCount;

if (debug == 9) {
    fprintf(stdout, "l1 ~ n3:%3d%3d%3d%3d%3d%3d\n", l1, m1, n1, l3, m3, n3);
    PrintBasis(b1, inp->gp, inp->gXYZ);
    PrintBasis(b2, inp->gp, inp->gXYZ);
    PrintBasis(b3, inp->gp, inp->gXYZ);
    PrintBasis(b4, inp->gp, inp->gXYZ);
}

    L = GetAngularMomentum(b1->Type) + GetAngularMomentum(b2->Type) + GetAngularMomentum(b3->Type) + GetAngularMomentum(b4->Type);

    F = calloc(sizeof(double), L + 1);

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            g1 = &b1->gaussian[i];
            g2 = &b2->gaussian[j];
            gid1 = g1->gid;
            gid2 = g2->gid;
            cid1 = inp->gp[gid1].cid;
            cid2 = inp->gp[gid2].cid;

            zeta = inp->zeta[gid1][gid2];

            PP(inp->P[gid1][gid2], inp->gXYZ[cid1], PA);
            PP(inp->P[gid1][gid2], inp->gXYZ[cid2], PB);

            KAB = inp->K[gid1][gid2];

            for (k = 0; k < gaussCount_3; k++) {
                for (l = 0; l < gaussCount_4; l++) {
                    g3 = &b3->gaussian[k];
                    g4 = &b4->gaussian[l];
                    gid3 = g3->gid;
                    gid4 = g4->gid;
                    cid3 = inp->gp[gid3].cid;
                    cid4 = inp->gp[gid4].cid;

                    eta = inp->zeta[gid3][gid4];

                    PP(inp->P[gid3][gid4], inp->gXYZ[cid3], QC);
                    PP(inp->P[gid3][gid4], inp->gXYZ[cid4], QD);
                    PP(inp->P[gid1][gid2], inp->P[gid3][gid4], PQ);

                    MULTISCALE(WP, PQ, -eta / (zeta + eta));
                    MULTISCALE(WQ, PQ, zeta / (zeta + eta));

                    rho = zeta * eta / (zeta + eta);

                    T = rho * NORM_2(PQ);

                    if (T < 17) {
                        for (m = 0; m <= L; m++)
                            F[m] = F_inc_gamma(m, T);
                    }else{
                        F[L] = F_inc_gamma(L, T);
                        double expT = exp(-T);
                        for (m = L - 1; m >= 0; m--)
                            F[m] = (2*T * F[m+1] + expT) / (2*m+1);
                    }
                    //fprintf(stdout, "%16.9E%16.9E%16.9E%16.9E%16.9E\n", rho, PQ.x, PQ.y, PQ.z, T);

                    KCD = inp->K[gid3][gid4];

                    pre = g1->coeff * g2->coeff * g3->coeff * g4->coeff;
/*
                    fprintf(stdout, "----- ----- ----- -----\n \
                                     %d %d %d\n \
                                     %d %d %d\n \
                                     %d %d %d\n \
                                     %d %d %d\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E %15.9E %15.9E\n \
                                     %15.9E\n", 
                                    g1->l, g1->m ,g1->n,
                                    g2->l, g2->m ,g2->n,
                                    g3->l, g3->m ,g3->n,
                                    g4->l, g4->m ,g4->n,
                                    zeta, eta, rho,
                                    pre, KAB, KCD,
                                    inp->P[gid1][gid2].x, inp->P[gid1][gid2].y, inp->P[gid1][gid2].z,
                                    inp->P[gid3][gid4].x, inp->P[gid3][gid4].y, inp->P[gid3][gid4].z,
                                    PA.x, PA.y, PA.z,
                                    PB.x, PB.y, PB.z,
                                    QC.x, QC.y, QC.z,
                                    QD.x, QD.y, QD.z,
                                    WQ.x, WQ.y, WQ.z,
                                    WP.x, WP.y, WP.z,
                                    T);
*/                                    
/* PASS */
/*
fprintf(stdout, "===== ===== =====\n%15.9E %15.9E %15.9E %15.9E\n",
                    inp->gp[gid1].alpha, inp->gp[gid2].alpha,
                    inp->gp[gid3].alpha, inp->gp[gid4].alpha);
                    fprintf(stdout, "----- ----- -----\n%15.9E %15.9E\n", zeta, eta);
*/

                    result += pre * KAB * KCD * HGPHrrVRR(l1, m1, n1, l3, m3, n3,
                                                       zeta, eta, rho,
                                                       &PA, &PB, &QC, &QD, &WQ, &WP,
                                                       0, F)/sqrt(zeta + eta);
                }
            }
        }
    }
    free(F);
    if (debug == 999)
        fprintf(stdout, "result: %15.8lf\n", result);
    return result;
}

double HGPHrrVRR(int l1, int m1, int n1, int l3, int m3, int n3,
            double zeta, double gamma, double ro,
            const COORD *PA, const COORD *PB, const COORD *QC,
            const COORD *QD, const COORD *WP, const COORD *WQ,
            int m, double *T)
{
// VRR
// compute the primitive integrals [e0|f0]
    double item1 = 0, item2 = 0, item3 = 0;
    double result;
    if (n3 >= 1) {
        if (fabs(QC->z) > COORDINATION_THRESHOLD)
            item1 += QC->z * HGPHrrVRR(l1, m1, n1, l3, m3, n3-1,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);
        if (fabs(WQ->z) > COORDINATION_THRESHOLD)
            item1 += WQ->z * HGPHrrVRR(l1, m1, n1, l3, m3, n3-1,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n3 >= 2) {
            item2 = (n3-1)/(2*gamma) * (HGPHrrVRR(l1, m1, n1, l3, m3, n3-2,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro/gamma * HGPHrrVRR(l1, m1, n1, l3, m3, n3-2,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n1 >= 1)
            item3 = n1/(2*(zeta + gamma))*HGPHrrVRR(l1, m1, n1-1, l3, m3, n3-1,
                        zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (m3 >= 1) {
        if (fabs(QC->y) > COORDINATION_THRESHOLD)
            item1 += QC->y * HGPHrrVRR(l1, m1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);
        if (fabs(WQ->y) > COORDINATION_THRESHOLD)
            item1 += WQ->y * HGPHrrVRR(l1, m1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        if (m3 >= 2) {
            item2 = (m3-1)/(2*gamma) * (HGPHrrVRR(l1, m1, n1, l3, m3-2, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * HGPHrrVRR(l1, m1, n1, l3, m3-2, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m1 >= 1)
            item3 = m1/(2*(zeta + gamma))*HGPHrrVRR(l1, m1-1, n1, l3, m3-1, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (l3 >= 1) {
        if (fabs(QC->x) > COORDINATION_THRESHOLD)
            item1 += QC->x * HGPHrrVRR(l1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);
        if (fabs(WQ->x) > COORDINATION_THRESHOLD)
            item1 += WQ->x * HGPHrrVRR(l1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        if (l3 >= 2) {
            item2 = (l3-1) / (2*gamma) * (HGPHrrVRR(l1, m1, n1, l3-2, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * HGPHrrVRR(l1, m1, n1, l3-2, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l1 >= 1)
            item3 = l1/(2*(zeta + gamma))*HGPHrrVRR(l1-1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }

    if (n1 >= 1) {
        if (fabs(PA->z) > COORDINATION_THRESHOLD)
            item1 += PA->z * HGPHrrVRR(l1, m1, n1-1, l3, m3, n3,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);
        if (fabs(WP->z) > COORDINATION_THRESHOLD)
            item1 += WP->z * HGPHrrVRR(l1, m1, n1-1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        if (n1 >= 2) {
            item2 = (n1-1) / (2*zeta) * (HGPHrrVRR(l1, m1, n1-2, l3, m3, n3,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * HGPHrrVRR(l1, m1, n1-2, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n3 >= 1)
            item3 = n3/(2*(zeta + gamma))*HGPHrrVRR(l1, m1, n1-1, l3, m3, n3-1,
                        zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (m1 >= 1) {
        if (fabs(PA->y) > COORDINATION_THRESHOLD)
            item1 += PA->y * HGPHrrVRR(l1, m1-1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);
        if (fabs(WP->y) > COORDINATION_THRESHOLD)
            item1 += WP->y * HGPHrrVRR(l1, m1-1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m1 >= 2) {
            item2 = (m1-1) / (2*zeta) * (HGPHrrVRR(l1, m1-2, n1, l3, m3, n3,
                                zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * HGPHrrVRR(l1, m1-2, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m3 >= 1)
            item3 = m3/(2*(zeta + gamma))*HGPHrrVRR(l1, m1-1, n1, l3, m3-1, n3,
                        zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }
    if (l1 >= 1) {
        if (fabs(PA->x) < COORDINATION_THRESHOLD)
            item1 += PA->x * HGPHrrVRR(l1-1, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);
        if (fabs(WP->x) < COORDINATION_THRESHOLD)
            item1 += WP->x * HGPHrrVRR(l1-1, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        if (l1 >= 2) {
            item2 = (l1-1) / (2*zeta) * (HGPHrrVRR(l1-2, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / zeta * HGPHrrVRR(l1-2, m1, n1, l3, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l3 > 0)
            item3 = l3/(2*(zeta + gamma))*HGPHrrVRR(l1-1, m1, n1, l3-1, m3, n3,
                            zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item3 = 0;

        result = item1 + item2 + item3;
        return result;
    }

    return T[m];
}
