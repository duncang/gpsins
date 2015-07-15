/**
 * \file matrix.h  Matrix routines for GPS and INS
 * $Id: matrix.h 955 2007-12-06 06:38:04Z greerd $
 */

#ifndef _MATRIX_H
#define _MATRIX_H

void MatMul331(const double Mat33[3][3], const double *Mat31, double *Result);
void MatMul333(const double Mat33a[3][3], const double Mat33b[3][3], double *Result);

double pythag (double a, double b);
/* singular value decomposition algorithm */
int svdcmp( int m, int n, double a[m][n], double w[n], double v[n][n], double work[n]);
void printvec(int n, double vec[n]);
void printm(int m,int n,double mat[m][n]);

#endif
