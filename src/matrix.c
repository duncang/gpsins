/**
 * \file matrix.c  Matrix routines for GPS and INS
 * $Id: matrix.c 955 2007-12-06 06:38:04Z greerd $
 */


#include <math.h>
#include <string.h> /* needed for memcpy */
#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))
static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
(dminarg1) : (dminarg2))
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))
static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
(minarg1) : (minarg2))
static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
(lmaxarg1) : (lmaxarg2))
static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
(lminarg1) : (lminarg2))
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
(imaxarg1) : (imaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))


void MatMul331(const double Mat33[3][3], const double *Mat31, double *Result)
{
	double LocalResult[3];
	
	LocalResult[0] = Mat33[0][0] * Mat31[0] + Mat33[0][1] * Mat31[1] + Mat33[0][2] * Mat31[2];
	LocalResult[1] = Mat33[1][0] * Mat31[0] + Mat33[1][1] * Mat31[1] + Mat33[1][2] * Mat31[2];
	LocalResult[2] = Mat33[2][0] * Mat31[0] + Mat33[2][1] * Mat31[1] + Mat33[2][2] * Mat31[2];
	
	memcpy(Result,LocalResult,sizeof(double)*3);
}

void MatMul333(const double Mat33a[3][3], const double Mat33b[3][3], double *Result)
{
	double LocalResult[3][3];
	int iRowIndex,iColIndex;
	
	for(iRowIndex = 0;iRowIndex<3;iRowIndex++)
	{
	  for(iColIndex=0;iColIndex<3;iColIndex++)
	  {
	    LocalResult[iRowIndex][iColIndex] =  Mat33a[iRowIndex][0]*Mat33b[0][iColIndex] +
						 Mat33a[iRowIndex][1]*Mat33b[1][iColIndex] +
						 Mat33a[iRowIndex][2]*Mat33b[2][iColIndex]; 
	  }
	}	
	memcpy(Result,LocalResult,sizeof(LocalResult));

}

/** 
	Computes (a2 + b2)1/2 without destructive underflow or overflow.
*/
double pythag (double a, double b)
{
  double absa, absb;
  absa = fabs (a);
  absb = fabs (b);
  if (absa > absb)
    return absa * sqrt (1.0 + DSQR(absb / absa));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt (1.0 + DSQR (absa / absb)));
}


/**
 * \fn svdcmp singular value decomposition algorithm 
 * \param a input matrix m by n
 * \param m number of rows
 * \param n number of columns
 * \param w vector length n containing diagonal elements of W
 * \param v n by n matrix V
 * \param work workspace vector of length n 
 *
 * \return -1 if iteration limit exceeded
 * 
 * The function calculates the singular value decomposition of matrix
 * A such that A = U * W * V'.  U replaces A on the output.  U is nominally 
 * a m by m matrix however if m < n, only the first m rows affect the 
 * result.  
 * 
 */
int svdcmp(const int m, const int n, double a[m][n], double w[n], double v[n][n], double work[n])
{


  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;

  int iReturnValue = 0;



  g=scale=anorm=0.0;  
  for (i=1;i<=n;i++) {
    l=i+1;
    work[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) 
        scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<=m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=i;k<=m;k++) 
            s += a[k][i]*a[k][j];
          f=s/h;
          for (k=i;k<=m;k++) 
            a[k][j] += f*a[k][i];
        }
        for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<=n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<=n;k++) work[k]=a[i][k]/h;
        for (j=l;j<=m;j++) {
          for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
          for (k=l;k<=n;k++) a[j][k] += s*work[k];
        }
        for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(fabs(w[i])+fabs(work[i])));
  }
  for (i=n;i>=1;i--) {    
    if (i < n) {
      if (g) {
        for (j=l;j<=n;j++)   
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=work[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) {    
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
        for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
        f=(s/a[i][i])*g;
        for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {    
    for (its=1;its<=30;its++) {   
      flag=1;
      for (l=k;l>=1;l--) {    
        nm=l-1;    
        if ((fabs(work[l])+anorm) == anorm) {
          flag=0;
          break;
        }
        if ((fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
        c=0.0;    
        s=1.0;
        for (i=l;i<=k;i++) {
          f=s*work[i];
          work[i]=c*work[i];
          if ((fabs(f)+anorm) == anorm) break;
          g=w[i];
          h=pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=1;j<=m;j++) {
            y=a[j][nm];
            z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l == k) {   
        if (z < 0.0) {   
          w[k] = -z;
          for (j=1;j<=n;j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 30) 
      { 
	iReturnValue = -1; 
      }
      x=w[l];    
      nm=k-1;
      y=w[nm];
      g=work[nm];
      h=work[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;    
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=work[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=pythag(f,h);
        work[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g = g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=1;jj<=n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j]=z;    
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=1;jj<=m;jj++) {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      work[l]=0.0;
      work[k]=f;
      w[k]=x;
    }
  }

  return iReturnValue;

}

/**
 Prints an mxn matrix
  note that this declaration will only work in gnu C..
   see http://gcc.gnu.org/onlinedocs/gcc-2.95.3/gcc_4.html#SEC75
 */

void printm(int m,int n,double mat[m][n])
{
    int i,j;
    for(i = 0; i < m; i++)
    {
    	for(j = 0; j < n; j++)
	{
	    fprintf(stdout, "%7.5f ", mat[i][j]);
	}
	fprintf(stdout,"\n");
    }			
    fprintf(stdout, "\n");
}

/* Prints a vector nx1 */
void printvec(int n, double vec[n])
{
    int i;
    for(i = 0; i < n; i++)
    {
	fprintf(stdout, "%7.5f ", vec[i]);
    }			
    fprintf(stdout, "\n");

}

 


