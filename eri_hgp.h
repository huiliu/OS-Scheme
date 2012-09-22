#include "common.h"


#ifndef ERI__HGP
#define ERI__HGP

double ERI_HGP_DRIVE(const INPUT_INFO *inp, int ii, int jj, int kk, int ll, int debug);
double HGPBasisHRR(const INPUT_INFO *inp,
                    int ii, int jj, int kk, int ll, 
                    int l1, int m1, int n1,
                    int l2, int m2, int n2,
                    int l3, int m3, int n3,
                    int l4, int m4, int n4,
                    COORD *AB, COORD *CD, double *XSXS, int debug);
double HGPBasis(const INPUT_INFO *inp,
                int ii, int jj, int kk, int ll,
                int l1, int m1, int n1,
                int l3, int m3, int n3,
                double *XSXS, int debug);
double HGPHrrVRR(int l1, int m1, int n1, int l3, int m3, int n3,
            double zeta, double gamma, double ro,
            const COORD *PA, const COORD *PB, const COORD *QC,
            const COORD *QD, const COORD *WP, const COORD *WQ,
            int m, double *T);

#endif
