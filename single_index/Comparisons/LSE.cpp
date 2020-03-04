//
//  ComputeLSE
//  Profile LS estimator for the single index model
//
//  Created by Piet Groeneboom on 02-02-20.
//  Copyright (c) 2018 Piet. All rights reserved.
//

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
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

int         m,n,nIterations;
double      **xx,*yy,*yy1,*vv,*f,*cumw,*cs;
double      *psi,*derivative,beta1,*alpha;

double  criterion(double beta1);
void    sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
int     CompareTime(const void *a, const void *b);
double  golden(double a1, double b1, double (*f)(double));
void    swap(double *x,double *y);


// [[Rcpp::export]]

List ComputeLSE(NumericMatrix X, NumericVector y, int m1)
{
    int     i,j;
    double  pi;
    
    pi=4*atan(1);
    
    beta1=pi/4;
       
    // determine the sample size
    
    n = (int)(y.size());
    
    // m is the dimension
    m= (int)m1;
    
    nIterations=100;
    
    // copy the data vector for use of the C++ procedures
    
    yy = new double[n];
    yy1 = new double[n+1];
    
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
    f= new double[m];
    
    cumw = new double[n+1];
    cs = new double[n+1];
    psi  = new double[n];
    
    cumw[0]=cs[0]=0;
    
    beta1 = golden(0,pi/2,criterion);
    
    alpha[0]=cos(beta1);
    alpha[1]=sin(beta1);

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
    
    delete[] yy, delete[] yy1; delete[] vv; delete[] alpha; delete[] f;
    delete[] cumw; delete[] cs; delete[] psi;
    
    return out;
}


double criterion(double beta1)
{
    int i;
    double sum;
    
    alpha[0]=cos(beta1);
    alpha[1]=sin(beta1);
    
    sort_alpha(m,n,xx,alpha,vv,yy);
    
    cs[0]=cumw[0]=0;
      
    yy1[0]=0;
      
    for (i=1;i<=n;i++)
        yy1[i]=yy[i-1];
      
    for (i=1;i<=n;i++)
    {
        cumw[i]=i*1.0;
        cs[i]=cs[i-1]+yy1[i];
    }
    
    convexmin(n,cumw,cs,yy1);
    
    for (i=0;i<n;i++)
        psi[i]=yy1[i+1];
    
    sum=0;
    
    for (i=0;i<n;i++)
        sum += SQR(psi[i]-yy[i])/n;
    
    return sum;
    
}

double golden(double a1, double b1, double (*f)(double))
{
  double a,b,eps=1.0e-10;
  
  a=a1;
  b=b1;
  
  double k = (sqrt(5.0) - 1.0) / 2;
  double xL = b - k*(b - a);
  double xR = a + k*(b - a);
  
  while (b-a>eps)
  {
    if ((*f)(xL)<(*f)(xR))
    {
      b = xR;
      xR = xL;
      xL = b - k*(b - a);
    }
    else
    {
      a = xL;
      xL = xR;
      xR = a + k * (b - a);
    }
  }
  return (a+b)/2;
  
}


void convexmin(int n, double cumw[], double cs[], double y[])
{
  int    i, j, m;

  y[1] = cs[1]/cumw[1];
  for (i=2;i<=n;i++)
  {
    y[i] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
    if (y[i-1]>y[i])
    {
      j = i;
      while (y[j-1] > y[i] && j>1)
      {
        j--;
        if (j>1)
          y[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
        else
          y[i] = cs[i]/cumw[i];
        for (m=j;m<i;m++)    y[m] = y[i];
      }
    }
  }
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


