#include <stdio.h>
#include <math.h>
#include "eri_os.h"


double ERI_basis_OS(const BASIS* b1, const BASIS* b2,
                    const BASIS* b3, const BASIS* b4, const INPUT_INFO *inp,int debug)
{ 
    int i, j, k, l, gaussCount_1, gaussCount_2, gaussCount_3, gaussCount_4;
    int gid1, gid2, gid3, gid4;
    int cid1, cid2, cid3, cid4;
    int L, m;
    GTO *g1, *g2, *g3, *g4;
    double zeta, eta, rho;
    COORD PA, PB, QC, QD, WQ, WP, PQ; 
    double KAB, KCD, pre;
    double *F;
    double T;
    double result = 0;

    gaussCount_1 = b1->gaussCount;
    gaussCount_2 = b2->gaussCount;
    gaussCount_3 = b3->gaussCount;
    gaussCount_4 = b4->gaussCount;
/*
    PrintBasis(b1, inp->gp, inp->gXYZ);
    PrintBasis(b2, inp->gp, inp->gXYZ);
    PrintBasis(b3, inp->gp, inp->gXYZ);
    PrintBasis(b4, inp->gp, inp->gXYZ);
*/

    L = GetAngularMomentum(b1->Type) + GetAngularMomentum(b2->Type) + GetAngularMomentum(b3->Type) + GetAngularMomentum(b4->Type);

    F = calloc(sizeof(double), L + 1);

    for (i = 0; i < gaussCount_1; i++) {
        for (j = 0; j < gaussCount_2; j++) {
            for (k = 0; k < gaussCount_3; k++) {
                for (l = 0; l < gaussCount_4; l++) {

                    g1 = &b1->gaussian[i];
                    g2 = &b2->gaussian[j];
                    g3 = &b3->gaussian[k];
                    g4 = &b4->gaussian[l];

                    gid1 = g1->gid;
                    gid2 = g2->gid;
                    gid3 = g3->gid;
                    gid4 = g4->gid;

                    cid1 = inp->gp[gid1].cid;
                    cid2 = inp->gp[gid2].cid;
                    cid3 = inp->gp[gid3].cid;
                    cid4 = inp->gp[gid4].cid;

                    zeta = inp->zeta[gid1][gid2];
                    eta = inp->zeta[gid3][gid4];

                    PP(inp->P[gid1][gid2], inp->gXYZ[cid1], PA);
                    PP(inp->P[gid1][gid2], inp->gXYZ[cid2], PB);
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

                    KAB = inp->K[gid1][gid2];
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
/* PASS
fprintf(stdout, "===== ===== =====\n%15.9E %15.9E %15.9E %15.9E\n",
                    inp->gp[gid1].alpha, inp->gp[gid2].alpha,
                    inp->gp[gid3].alpha, inp->gp[gid4].alpha);
                    fprintf(stdout, "----- ----- -----\n%15.9E %15.9E\n", zeta, eta);
*/

                    result += pre * KAB * KCD * ERI_VRR_OS(g1->l, g1->m ,g1->n,
                                                       g2->l, g2->m ,g2->n,
                                                       g3->l, g3->m ,g3->n,
                                                       g4->l, g4->m ,g4->n,
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

#define COORDINATION_THRESHOLD          1.0E-12
double ERI_VRR_OS(int l1, int m1, int n1,
                  int l2, int m2, int n2,
                  int l3, int m3, int n3,
                  int l4, int m4, int n4,
                  double zeta, double gamma, double ro,
                  const COORD *PA, const COORD *PB, const COORD *QC,
                  const COORD *QD, const COORD *WQ, const COORD *WP,
                  int m, double *T)
{
    double item1, item11, item12, item2, item3, item4, item5, result;
    // -----------------------------------the fourth GTO-------------------------------
    if (n4 >= 1) {
/*
        item1 = QD->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(QD->z) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = QD->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WQ->z) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WQ->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n4 >= 2) {
            item2 = (n4-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-2, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4-2, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n3 >=1)
            item3 = n3 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n2 >= 1)
            item4 = n2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (n1 >= 1)
            item5 = n1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m4 >= 1) {
/*
        item1 = QD->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(QD->y) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = QD->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WQ->y) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WQ->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m4 >= 2) {
            item2 = (m4-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-2, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4-2, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m3 >=1)
            item3 = m3 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m2 >= 1)
            item4 = m2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (m1 >= 1)
            item5 = m1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l4 >= 1) {
/*
        item1 = QD->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(QD->x) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = QD->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WQ->x) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WQ->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l4 >= 2) {
            item2 = (l4-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-2, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3, l4-2, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l3 >=1)
            item3 = l3 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T)
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l2 >= 1)
            item4 = l2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (l1 >= 1)
            item5 = l1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }
    // -----------------------------------the third GTO-------------------------------
    if (n3 >= 1) {
/*
        item1 = QC->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(QC->z) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = QC->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WQ->z) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WQ->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n3 >= 2) {
            item2 = (n3-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-2, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-2, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n4 >=1)
            item3 = n4 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n1 >= 1)
            item4 = n1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (n2 >= 1)
            item5 = n2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m3 >= 1) {
/*
        item1 = QC->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(QC->y) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = QC->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WQ->y) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WQ->y * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m3 >= 2) {
            item2 = (m3-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-2, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-2, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m4 >=1)
            item3 = m4 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m1 >= 1)
            item5 = m1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (m2 >= 1)
            item4 = m2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l3 >= 1) {
/*
        item1 = QC->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WQ->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(QC->x) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = QC->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WQ->x) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WQ->x * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l3 >= 2) {
            item2 = (l3-1) / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-2, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-2, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l4 >=1)
            item3 = l4 / (2*gamma) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / gamma * ERI_VRR_OS(l1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l1 >= 1)
            item5 = l1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (l2 >= 1)
            item4 = l2 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }
    // ---------------------------------the second basis-------------------------------------
    if (n2 >= 1) {
/*
        item1 = PB->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(PB->z) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = PB->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WP->z) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WP->z * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n2 >= 2) {
            item2 = (n2-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1, l2, m2, n2-2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n1 >=1)
            item3 = n1 / (2*zeta) * (ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n4 >= 1)
            item4 = n4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP,  m+1, T);
        else
            item4 = 0;

        if (n3 >= 1)
            item5 = n1 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2, n2-1, l3, m3, n3-1, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP,  m+1, T);
        else
            item5 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m2 >= 1) {
/*
        item1 = PB->y * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->y * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(PB->y) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = PB->y * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WP->y) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WP->y * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m2 >= 2) {
            item2 = (m2-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1, l2, m2-2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1, l2, m2-2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m1 >=1)
            item3 = m1 / (2*zeta) * (ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m4 >= 1)
            item5 = m4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (m3 >= 1)
            item4 = m3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2, m2-1, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l2 >= 1) {
/*
        item1 = PB->x * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->x * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(PB->x) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = PB->x * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WP->x) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WP->x * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l2 >= 2) {
            item2 = (l2-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1, l2-2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1, l2-2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l1 >=1)
            item3 = l1 / (2*zeta) * (ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l4 >= 1)
            item5 = l4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (l3 >= 1)
            item4 = l3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1, l2-1, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    // the first basis
    if (n1 >= 1) {
/*
        item1 = PA->z * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->z * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(PA->z) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = PA->z * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WP->z) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WP->z * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (n1 >= 2) {
            item2 = (n1-1) / (2*zeta) * (ERI_VRR_OS(l1, m1, n1-2, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1-2, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (n2 >=1)
            item3 = n2 / (2*zeta) * (ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2-1, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (n3 >= 1)
            item4 = n3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3-1, l4, m4, n4,  zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        if (n4 >= 1)
            item5 = n4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1, n1-1, l2, m2, n2, l3, m3, n3, l4, m4, n4-1, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP,  m+1, T);
        else
            item5 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (m1 >= 1) {
/*
        item1 = PA->y * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->y * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/
        if (fabs(PA->y) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = PA->y * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WP->y) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WP->y * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (m1 >= 2) {
            item2 = (m1-1) / (2*zeta) * (ERI_VRR_OS(l1, m1-2, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1-2, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (m2 >=1)
            item3 = m2 / (2*zeta) * (ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1, m1-1, n1, l2, m2-1, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (m3 >= 1)
            item5 = m3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3-1, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (m4 >= 1)
            item4 = m4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1, m1-1, n1, l2, m2, n2, l3, m3, n3, l4, m4-1, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    if (l1 >= 1) {
/*
        item1 = PA->x * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
              + WP->x * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
*/        
        if (fabs(PA->x) < COORDINATION_THRESHOLD)
            item11 = 0;
        else
            item11 = PA->x * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T);

        if (fabs(WP->x) < COORDINATION_THRESHOLD)
            item12 = 0;
        else
            item12 = WP->x * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);

        if (l1 >= 2) {
            item2 = (l1-1) / (2*zeta) * (ERI_VRR_OS(l1-2, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1-2, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        }else
            item2 = 0;

        if (l2 >=1)
            item3 = l2 / (2*zeta) * (ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m, T) \
                       - ro / zeta * ERI_VRR_OS(l1-1, m1, n1, l2-1, m2, n2, l3, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T));
        else
            item3 = 0;

        if (l3 >= 1)
            item5 = l3 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3-1, m3, n3, l4, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item5 = 0;

        if (l4 >= 1)
            item4 = l4 / (2*(zeta + gamma)) * ERI_VRR_OS(l1-1, m1, n1, l2, m2, n2, l3, m3, n3, l4-1, m4, n4, zeta, gamma, ro, PA, PB, QC, QD, WQ, WP, m+1, T);
        else
            item4 = 0;

        //result = item1 + item2 + item3 + item4 + item5;
        result = item11 + item12 + item2 + item3 + item4 + item5;
        return result;
    }

    return T[m];
}
