//
//  Compute spline
//  Spline estimation for the single index model
//
//  Created by Piet Groeneboom on 04-05-18.
//  Copyright (c) 2018 Piet Groeneboom. All rights reserved.
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

int m,n;
double **xx,*yy,*vv,*f,*psi,mu,mu1,**rr,penalty;
double **q,*d,*D,**L,**b,**b_inv,*derivative,*alpha,beta1;


double  criterion(double beta);
void    sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[]);
int     CompareTime(const void *a, const void *b);
void    swap(double *x,double *y);
void    Compute_Q(int n, double q[], double b[]);
void    Compute_cubic_spline();
void    Cholesky_sol(int n, double b[], double z[], double x[]);
void    Cholesky_dec(int n, double b[], double D[], double **L);
double  golden(double a1, double b1, double (*f)(double));
double  Compute_penalty(double gamma[]);
void    Compute_R();

// [[Rcpp::export]]

List Compute_spline(NumericMatrix X, NumericVector y, int m1)
{
    int  i,j;
    double pi;
    
    pi=4*atan(1);
       
    // determine the sample size
    
    n = (int)(y.size());
    
    // m is the dimension
    m= (int)m1;
    mu=0.1;
    
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
    f= new double[m];
    
    psi  = new double[n];
    derivative  = new double[n];
    
    rr = new double *[n];
    for (i=0;i<n;i++)
          rr[i] = new double [2];
      
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
    
    NumericMatrix out3 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out3(i,0)=vv[i];
        out3(i,1)=derivative[i];
    }
        
    // make the list for the output, containing alpha and the estimate of psi
        
    List out = List::create(Rcpp::Named("data")=out0,Rcpp::Named("alpha")=out1,Rcpp::Named("psi")=out2,
                                        Rcpp::Named("derivative")=out3);
    // free memory
   
    for (i=0;i<n;i++)
        delete[] xx[i];
    
    delete[] xx;
    
    for (i=0;i<n;i++)
           delete[] rr[i];
       
    delete[] rr;
    
    delete[] yy; delete[] vv; delete[] alpha; delete[] f; delete[] psi;
    delete[] derivative;
    
    return out;
}

// This is the criterion with the derivative of the spline

double criterion(double beta)
{
    int i;
    double sum;
    
    alpha[0]=cos(beta);
    alpha[1]=sin(beta);
    
    sort_alpha(m,n,xx,alpha,vv,yy);
    
    Compute_cubic_spline();
    
    sum=0;
    
    for (i=0;i<n;i++)
        sum += SQR(psi[i]-yy[i])/n;
    
    return sum+mu*penalty;
}

void Compute_R()
{
    int i;
    
    for (i=1;i<=n-2;i++)
    {
        rr[i][0]=(vv[i+1]-vv[i-1])/3;
        if (i<n-2)
            rr[i][1]=(vv[i+1]-vv[i])/6;
    }
}

double Compute_penalty(double gamma[])
{
    int i;
    double sum;
    
    Compute_R();
    
    sum=0;
    for (i=1;i<=n-2;i++)
    {
        sum += rr[i][0]*SQR(gamma[i]);
        if (i<n-2)
            sum += 2*rr[i][1]*gamma[i]*gamma[i+1];
    }
    
    return sum/n;
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

void Compute_cubic_spline()
{
    int i,j,m;
    double *a,*b,*c,*d,*q,*gamma;
    
    m=3*n;
    
    a= new double[m+1];
    b= new double[m+1];
    c= new double[n+1];
    d= new double[n+1];
    
    q= new double[m+1];
    
    gamma= new double[n];
    
    for (i=0;i<=m;i++)
        a[i]=b[i]=0;
    
    for (i=0;i<=n;i++)
        c[i]=d[i]=0;
    
    
    for (i=0;i<=m;i++)
        q[i]=0;
    
    
    Compute_Q(n,q,b);
    
    for (j=2;j<=n-1;j++)
        a[2*(j-2)+j-1] += (vv[j]-vv[j-2])/3;
    
    for (j=2;j<n-1;j++)
        a[2*(j-1)+j-1] += (vv[j]-vv[j-1])/6;
    
    
    for (i=1;i<=m;i++)
        a[i] += mu*b[i];
    
    for (i=2;i<=n-1;i++)
        c[i-1]=(yy[i]-yy[i-1])/(vv[i]-vv[i-1])-(yy[i-1]-yy[i-2])/(vv[i-1]-vv[i-2]);
    
    Cholesky_sol(n-2,a,c,gamma);
    
    penalty = Compute_penalty(gamma);
    
    d[1] += q[1]*gamma[1];
    
    d[2] += q[2]*gamma[1]+q[3]*gamma[2];
    
    for (i=3;i<=n-2;i++)
    {
        for (j=i-2;j<=i;j++)
            d[i] += q[i+2*(j-1)]*gamma[j];
    }
    
    for (j=n-3;j<=n-2;j++)
        d[n-1] += q[n-1+2*(j-1)]*gamma[j];
    
    d[n] += q[n+2*(n-3)]*gamma[n-2];
    
    for (i=0;i<n;i++)
        psi[i]=yy[i]-mu*d[i+1];
    
    derivative[0]=(psi[1]-psi[0])/(vv[1]-vv[0])-(1.0/6)*(vv[1]-vv[0])*gamma[1];
    derivative[n-1]=(psi[n-1]-psi[n-2])/(vv[n-1]-vv[n-2])+(1.0/6)*(vv[n-1]-vv[n-2])*gamma[n-2];
    
    for (j=1;j<=n-2;j++)
        derivative[j]=(psi[j]-psi[j-1])/(vv[j]-vv[j-1])+(1.0/6)*(vv[j]-vv[j-1])*(gamma[j-1]+2*gamma[j]);
    
    delete[] a; delete[] b; delete[] c; delete[] d;
    delete[] q; delete[] gamma;
    
}


void Compute_Q(int n, double q[], double b[])
{
    int i,j,k;
    
    for (j=2;j<=n-1;j++)
    {
        q[j-1+2*(j-2)] += 1.0/(vv[j-1]-vv[j-2]);
        q[j+2*(j-2)] += -1.0/(vv[j-1]-vv[j-2])-1.0/(vv[j]-vv[j-1]);
        q[j+1+2*(j-2)] += 1.0/(vv[j]-vv[j-1]);
    }
    
    
    for (i=2;i<=n-1;i++)
    {
        for (j=i;j<=i+2;j++)
        {
            for (k=j-1;k<=i+1;k++)
            {
                if (k>=1)
                    b[i-1+2*(j-2)] += q[k+2*(i-2)]*q[k+2*(j-2)];
            }
        }
    }
    
    
}


// The following two routine implement section 2.6.1 and 2.6.2 of Green and Silverman (1994)
//

void Cholesky_dec(int n, double b[], double D[], double **L)
{
    int i;
    
    // b is an an array version of the matrix B in (2.29) on p. 25 of Green and Silverman (1994)
    // B[i][j] = b[(2*(i-1)+j]
    
    D[1]=b[1];
    L[2][1] = b[3]/D[1];
    D[2]= b[4]-SQR(L[2][1])*D[1];
    
    for (i=3;i<=n;i++)
    {
        L[i][2] = b[2*(i-1)+i-2]/D[i-2];
        L[i][1] = (b[2*(i-1)+i-1]-L[i-1][1]*L[i][2]*D[i-2])/D[i-1];
        D[i]=b[2*(i-1)+i]-SQR(L[i][1])*D[i-1]-SQR(L[i][2])*D[i-2];
    }
    
}


void Cholesky_sol(int n, double b[], double z[], double x[])
{
    int i,j;
    double **L,*D;
    double *u,*v;
    
    
    u= new double[n+1];
    v= new double[n+1];
    
    D= new double[n+1];
    
    L = new double *[n+1];
    for (i=0;i<n+1;i++)
        L[i] = new double [3];
    
    for (i=1;i<=n;i++)
    {
        D[i]=0;
        for (j=0;j<=2;j++)
            L[i][j]=0;
    }
    
    Cholesky_dec(n,b,D,L);
    
    u[1]=z[1];
    u[2]=z[2]-L[2][1]*u[1];
    
    for (i=3;i<=n;i++)
        u[i]=z[i]-L[i][1]*u[i-1]-L[i][2]*u[i-2];
    
    for (i=1;i<=n;i++)
        v[i]=u[i]/D[i];
    
    x[n]=v[n];
    x[n-1]=v[n-1]-L[n][1]*x[n];
    
    for (i=n-2;i>=1;i--)
        x[i]=v[i]-L[i+1][1]*x[i+1]-L[i+2][2]*x[i+2];
    
    delete[] u;
    delete[] v;
    
    delete[] D;
    
    for (i=0;i<n+1;i++)
        delete[] L[i];
    
    delete[] L;
    
}

double best_nearby (int m, double delta[], double point[], double prevbest,
                    double f(int m, double alpha[]), int *funevals)
{
    double ftmp,minf,*z;
    int i;
    
    z = new double[m];
    
    minf = prevbest;
    
    for ( i = 0; i < m; i++ )
        z[i] = point[i];
    
    for ( i = 0; i < m; i++ )
    {
        z[i] = point[i] + delta[i];
        
        ftmp = f(m,z);
        *funevals = *funevals + 1;
        
        if ( ftmp < minf )
            minf = ftmp;
        else
        {
            delta[i] = - delta[i];
            z[i] = point[i] + delta[i];
            ftmp = f(m,z);
            *funevals = *funevals + 1;
            
            if ( ftmp < minf )
                minf = ftmp;
            else
                z[i] = point[i];
        }
    }
    
    for ( i = 0; i < m; i++ )
        point[i] = z[i];
    
    delete [] z;
    
    return minf;
}

int hooke(int m, double startpt[], double endpt[], double rho, double eps,
          int itermax, double f(int m, double alpha[]))
{
    double *delta,fbefore;
    int i,iters,keep,funevals,count;
    double newf,*newx,steplength,tmp;
    bool verbose = false;
    double *xbefore;
    
    delta = new double[m];
    newx = new double[m];
    xbefore = new double[m];
    
    for ( i = 0; i < m; i++ )
        xbefore[i] = newx[i] = startpt[i];
    
    for ( i = 0; i < m; i++ )
    {
        if ( startpt[i] == 0.0 )
            delta[i] = rho;
        else
            delta[i] = rho*fabs(startpt[i]);
    }
    
    funevals = 0;
    steplength = rho;
    iters = 0;
    
    
    fbefore = f(m,newx);
    funevals = funevals + 1;
    newf = fbefore;
    
    while ( iters < itermax && eps < steplength )
    {
        iters = iters + 1;
        
        if (verbose)
        {
            cout << "\n";
            cout << "  FUNEVALS, = " << funevals
            << "  F(X) = " << fbefore << "\n";
            
            for ( i = 0; i < m; i++ )
            {
                cout << "  " << i + 1
                << "  " << xbefore[i] << "\n";
            }
        }
        //
        //  Find best new alpha, one coordinate at a time.
        //
        for ( i = 0; i < m; i++ )
            newx[i] = xbefore[i];
        
        
        
        newf = best_nearby(m,delta,newx,fbefore,f,&funevals);
        //
        //  If we made some improvements, pursue that direction.
        //
        keep = 1;
        count=0;
        
        while (newf<fbefore && keep == 1 && count<=100)
        {
            count++;
            for ( i = 0; i < m; i++ )
            {
                //
                //  Arrange the sign of DELTA.
                //
                if ( newx[i] <= xbefore[i] )
                    delta[i] = - fabs(delta[i]);
                else
                    delta[i] = fabs(delta[i]);
                //
                //  Now, move further in this direction.
                //
                tmp = xbefore[i];
                xbefore[i] = newx[i];
                newx[i] = newx[i] + newx[i] - tmp;
            }
            
            fbefore = newf;
            
            newf = best_nearby(m,delta,newx,fbefore,f,&funevals);
            //
            //  If the further (optimistic) move was bad...
            //
            if (fbefore <= newf)
                break;
            //
            //  Make sure that the differences between the new and the old points
            //  are due to actual displacements; beware of roundoff errors that
            //  might cause NEWF < FBEFORE.
            //
            keep = 0;
            
            for ( i = 0; i < m; i++ )
            {
                if ( 0.5 * fabs(delta[i]) < fabs(newx[i]-xbefore[i]))
                {
                    keep = 1;
                    break;
                }
            }
        }
        
        if (eps <= steplength && fbefore <= newf)
        {
            steplength = steplength * rho;
            for ( i = 0; i < m; i++ )
                delta[i] = delta[i] * rho;
        }
        
    }
    
    for ( i = 0; i < m; i++ )
        endpt[i] = xbefore[i];
    
    delete [] delta;
    delete [] newx;
    delete [] xbefore;
    
    return iters;
}

void swap(double *x,double *y)
{
    double temp;
    temp=*x;
    *x=*y;
    *y=temp;
}


