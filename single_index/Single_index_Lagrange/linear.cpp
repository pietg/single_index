//
//  ComputeSSE
//  Simple score estimator for the single index model
//
//  Created by Piet Groeneboom on 04-05-18.
//  Copyright (c) 2018 Piet. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <random>
#include <Rcpp.h>

#define SQR(x) ((x)*(x))

using namespace std;
using namespace Rcpp;

#define SQR(x) ((x)*(x))

typedef struct
{
    int index;
    double v;
    double y;
}
data_object;

int m,n;
double **xx,*yy,*vv,*psi,*x_mean;

double  criterion(double alpha[]);
void    sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
int     CompareTime(const void *a, const void *b);
void    compute_alpha(double alpha[]);
int     LUPDecompose(double **A, int N, double Tol, int *P);
void    LUPSolve(double **A, int *P, double *b, int N, double *x);
void    LUPInvert(double **A, int *P, int N, double **IA);
double  LUPDeterminant(double **A, int *P, int N);
void    swap(double *x,double *y);


// [[Rcpp::export]]

List Compute_linear(NumericMatrix X, NumericVector y, NumericVector alpha0, int m1)
{
    int             i,j;
    double          sum,*alpha;
       
    // determine the sample size
    n = (int)(y.size());
    
    // m is the dimension
    m= (int)m1;
    
    // copy the data vector for use of the C++ procedures
    
    yy = new double[n];
    
    for (i=0;i<n;i++)
        yy[i]=(double)y[i];
    
    xx = new double *[n];
    for (i=0;i<n;i++)
        xx[i] = new double [m];
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<m;j++)
          xx[i][j]=(double)X(i,j);
    }
    
    vv= new double[n];
    alpha= new double[m];

    psi  = new double[n];
    x_mean = new double[m];
    
    compute_alpha(alpha);
        
    sum=0;
    
    for (i=0;i<m;i++)
        sum += SQR(alpha[i]);
    
    for (i=0;i<m;i++)
        alpha[i]/=sqrt(sum);
    
    sort_alpha(m,n,xx,alpha,vv,yy);
    
    NumericMatrix out0 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out0(i,0)=vv[i];
        out0(i,1)=yy[i];
    }
    
    NumericVector out1 = NumericVector(m);
    
    // computation of alpha
    
    for (i=0;i<m;i++)
        out1(i)=alpha[i];
    
    NumericMatrix out2 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out2(i,0)=vv[i];
        out2(i,1)=psi[i];
    }
    
    // make the list for the output, containing alpha and the estimate of psi
    
    List out = List::create(Rcpp::Named("data")=out0,Rcpp::Named("alpha")=out1,Rcpp::Named("psi")=out2);
    
    
    // free memory
   
    for (i=0;i<n;i++)
        delete[] xx[i];
    
    delete[] xx;
    
    delete[] yy, delete[] vv; delete[] alpha;
    delete[] psi; delete[] x_mean;
    
    return out;
}

void compute_alpha(double alpha[])
{
    int i,j,k,*P;
    double y_mean,*b,**S,*f,sum,tol=1.0e-10;
    
    S = new double *[m];
    for (i=0;i<m;i++)
        S[i] = new double [m];
    
    b = new double[m];
    P = new int[m+1];
    f = new double[m];
    
    for (i=0;i<m;i++)
    {
        for (j=0;j<m;j++)
            S[i][j]=0;
    }
    
    for (j=0;j<m;j++)
    {
        sum=0;
        for (i=0;i<n;i++)
            sum += xx[i][j];
        x_mean[j] = sum/n;
    }
    
    y_mean=0;
    for (i=0;i<n;i++)
        y_mean += yy[i];
    
    y_mean/=n;
    
    
    for (i=0;i<m;i++)
    {
        for (j=0;j<m;j++)
        {
            for (k=0;k<n;k++)
                S[i][j] += (xx[k][i]-x_mean[i])*(xx[k][j]-x_mean[j]);
            S[i][j]/=n;
        }
    }
    
    for (i=0;i<m;i++)
        b[i]=0;
    
    for (i=0;i<m;i++)
    {
        for (j=0;j<n;j++)
            b[i] += (xx[j][i]-x_mean[i])*yy[j];
        b[i]/=n;
    }
    
    LUPDecompose(S,m,tol,P);
    LUPSolve(S,P,b,m,alpha);
    
    sort_alpha(m,n,xx,alpha,vv,yy);
    
    for (i=0;i<n;i++)
    {
        psi[i]=0;
        for (j=0;j<m;j++)
            psi[i] += alpha[j]*(xx[i][j]-x_mean[j]);
        psi[i] += y_mean;
    }
    
    for (i=0;i<m;i++)
        delete[] S[i];
    delete[] S;
    
    delete[] b; delete[] P; delete[] f;
}



int LUPDecompose(double **A, int N, double Tol, int *P)
{
    int i, j, k, imax;
    double maxA, *ptr, absA;
    
    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    
    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;
        
        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }
        
        if (maxA < Tol) return 0; //failure, matrix is degenerate
        
        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            
            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;
            
            //counting pivots starting from N (for determinant)
            P[N]++;
        }
        
        for (j = i + 1; j < N; j++)
        {
            A[j][i] /= A[i][i];
            
            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }
    
    return 1;  //decomposition done
}

// INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
// OUTPUT: x - solution vector of A*x=b

void LUPSolve(double **A, int *P, double *b, int N, double *x)
{
    for (int i = 0; i < N; i++)
    {
        x[i] = b[P[i]];
        
        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }
    
    for (int i = N - 1; i >= 0; i--)
    {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];
        
        x[i] = x[i] / A[i][i];
    }
}

// INPUT: A,P filled in LUPDecompose; N - dimension
// OUTPUT: IA is the inverse of the initial matrix

void LUPInvert(double **A, int *P, int N, double **IA)
{
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            if (P[i] == j)
                IA[i][j] = 1.0;
            else
                IA[i][j] = 0.0;
            
            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }
        for (int i = N - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
            
            IA[i][j] = IA[i][j] / A[i][i];
        }
    }
    
}

// INPUT: A,P filled in LUPDecompose; N - dimension.
// OUTPUT: Function returns the determinant of the initial matrix

double LUPDeterminant(double **A, int *P, int N)
{
    
    double det = A[0][0];
    
    for (int i = 1; i < N; i++)
        det *= A[i][i];
    
    if ((P[N] - N) % 2 == 0)
        return det;
    else
        return -det;
}

int CompareTime(const void *a, const void *b)
{
    if ((*(data_object *) a).v < (*(data_object *) b).v)
        return -1;
    if ((*(data_object *) a).v > (*(data_object *) b).v)
        return 1;
    return 0;
}

void sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[])
{
    int i,j,*ind;
    double **xx_new;
    data_object *obs;
    
    obs= new data_object[n];
    ind= new int[n];
    
    xx_new = new double *[n];
    for (i=0;i<n;i++)
        xx_new[i] = new double [m];
    
    for (i=0;i<n;i++)
    {
        vv[i]=0;
        for (j=0;j<m;j++)
            vv[i] += alpha[j]*xx[i][j];
    }
    
    for (i=0;i<n;i++)
    {
        obs[i].index=i;
        obs[i].v=vv[i];
        obs[i].y=yy[i];
    }
    
    qsort(obs,n,sizeof(data_object),CompareTime);
    
    for (i=0;i<n;i++)
        ind[i]=obs[i].index;
    
    
    for (i=0;i<n;i++)
        for (j=0;j<m;j++)
            xx_new[i][j]=xx[ind[i]][j];
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<m;j++)
            xx[i][j]=xx_new[i][j];
        vv[i]=obs[i].v;
        yy[i]=obs[i].y;
    }
    
    delete[] obs;
    
    delete[] ind;
    for (i=0;i<n;i++)
        delete[] xx_new[i];
    delete[] xx_new;
}

void swap(double *x,double *y)
{
    double temp;
    temp=*x;
    *x=*y;
    *y=temp;
}


