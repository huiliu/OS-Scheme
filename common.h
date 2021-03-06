/*
 * File: parsebasis.h
 * --------------------------------
 * This Program parse the basis function information which come from ESML.
 * Author:  liuhui
 * Date: Sat Sep  8 12:51:42 CST 2012
 */

#include <stdlib.h>

#ifndef __INTEGRAL__COMMON__
#define __INTEGRAL__COMMON__

/* 
 *  BASIC DATA STRUCTURE
 * --------------------------------------------------------
 */

// angstrom to bohr radius
#define ANGS_2_BOHR     1.88971616463

typedef struct {
    double x, y, z;
}COORD;

// Store the exponent and coefficient
typedef struct {
    double alpha;
    double coeff;
    unsigned short cid;     // direct the coordination
}GTO_PARTIAL;

typedef struct gto {
    short int l, m, n;      // store the angular momentum
    double coeff;           // store coefficient * normalization
    unsigned short gid;     // direct the element of exponent and coef
}GTO;

typedef struct _b {
// the parameter of basis function
    unsigned short int gaussCount;  // the number of primitive function
    unsigned short Type;
    double  scale;
    GTO* gaussian;
    unsigned short cid;             // direct the coordination
}BASIS;

typedef struct atom_INFORMATON_ {
// 定义一个原子信息, 包括核电荷数，基函数数目，元素符号，坐标，基函数
    unsigned short int n;                   // Atomic Number
    unsigned short int basisCount;          // basis function count
    char symbol[3];                         // Element Symbol
    unsigned short cid;                     // direct the coordination
}ATOM_INFO;

typedef struct FILE_INPUT_ {
// 存储解析后的整个输入文件
    unsigned short int atomCount;           // atom count in current system
    unsigned short int eCount;              // electron count in the 
    short int icharge;                      // 
    unsigned short int imult;               // multi degree
    unsigned short int basisCount;          // basis function count
    unsigned short int gCount;            // independent GTO count
    unsigned short int gtoCount;            // independent GTO count
    BASIS *basisSet;                        // save all of the basisset
    GTO_PARTIAL *gp;
    GTO *gtoSet;
    COORD *gXYZ;                           // save all of coordination
    ATOM_INFO* atomList;                   // element information
    // --------------------------
    COORD **P;
    double **K;
    double **zeta;
}INPUT_INFO;


/*
 *  Function Imeplement
 * -----------------------------------------------------------------------
 */

inline double F_inc_gamma(int m ,double w);
inline int factorial(int);
inline int factorial_2(int);
// check the symtery of two-electron integral
inline int ChkERISym(double ****e, int i, int j, int k, int l, int N,
                                                               int *is_dup);
inline void Gaussian_product_center(const GTO_PARTIAL *gpa,
                            const GTO_PARTIAL *gpb, const COORD *c, COORD *p);

inline int GetAngularMomentum(int);
#define MALLOC(p,n) \
    if (!(p = malloc((n)))) { \
        printf("====== alloc memory failed! ======\n");\
        exit(EXIT_FAILURE);\
    }

#define CALLOC(p,n,s) \
    if (!(p = calloc(n, s))) { \
        printf("====== alloc memory failed! ======\n");\
        exit(EXIT_FAILURE);\
    }

#define REALLOC(p,n) \
    if (!(realloc(p, n))) { \
        printf("====== realloc memory failed! ======\n");\
        exit(EXIT_FAILURE);\
    }

#define DISTANCE(A,B)   (pow((A.x-B.x),2)+pow((A.y-B.y),2)+pow((A.z-B.z),2))
// The Norm of vector.
#define NORM_2(A)   (A.x*A.x + A.y*A.y + A.z*A.z)

#define MULTISCALE(A,B,a) \
    A.x = a * B.x; \
    A.y = a * B.y; \
    A.z = a * B.z

#define PP(A,B,P)   \
    P.x = A.x - B.x; \
    P.y = A.y - B.y; \
    P.z = A.z - B.z;

#endif
