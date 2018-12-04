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
double **xx,*yy,*vv,*cumw,*cs,*f,*psi;

double  criterion(double alpha[]);
void    sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
int     CompareTime(const void *a, const void *b);
void    nelmin ( double fn ( double x[] ), int n, double start[], double xmin[],
                double *ynewlo, double reqmin, double step[], int konvge, int kcount,
                int *icount, int *numres, int *ifault );
double  K(double x);
void    swap(double *x,double *y);


// [[Rcpp::export]]

List ComputeSSE(NumericMatrix X, NumericVector y, NumericVector alpha0, int m1)
{
    int             i,j,NumIt,seed;
    double          sum,*alpha,*alpha_init;
    double          *step,reqmin,ynewlo;
    int             icount,ifault,kcount,konvge,numres;
    
    // determine the sample size
    n = (int)(y.size());
    
    // m is the dimension
    m= (int)m1;
    
    seed=180;
    NumIt=20;
    
    std::mt19937 generator(seed);
    std::uniform_real_distribution<> dis(0,1);
    std::normal_distribution<double> dis_normal(0.0,1.0);
    
    reqmin=1.0e-8;
    konvge=10;
    kcount=2000;
    
    step = new double[m];
    
    for (i=0;i<m;i++)
        step[i]=1.0;
    
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
    alpha_init= new double[m];
    f= new double[m];
    
    for (i=0;i<m;i++)
        alpha_init[i]=(double)alpha0(i);
    
    cumw = new double[n+1];
    cs = new double[n+1];
    psi  = new double[n];
    
    cumw[0]=cs[0]=0;
    
    nelmin(criterion,m,alpha_init,alpha,&ynewlo,reqmin,step,konvge,kcount,&icount,&numres,&ifault);
        
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
    
    delete[] yy, delete[] vv; delete[] alpha; delete[] alpha_init; delete[] f;
    delete[] cumw; delete[] cs; delete[] psi; delete[] step ;
    
    return out;
}

double criterion(double alpha[])
{
    int i,j;
    double sum,lambda;

    sort_alpha(m,n,xx,alpha,vv,yy);
    
    for (i=1;i<=n;i++)
    {
        cumw[i]=i*1.0;
        cs[i]=cs[i-1]+yy[i-1];
    }
    
    convexmin(n,cumw,cs,psi);
    
    for (j=0;j<m;j++)
        f[j]=0;
    
    for (j=0;j<m;j++)
    {
        for (i=0;i<n;i++)
            f[j] += xx[i][j]*(psi[i]-yy[i]);
    }
    
    lambda=0;
    
    for (i=0;i<m;i++)
        lambda += alpha[i]*f[i];
    
    for (i=0;i<m;i++)
        f[i] -= lambda*alpha[i];
    
    sum=0;
    
    for (i=0;i<m;i++)
        sum += SQR(f[i]);
    
    return sum;
    
}


void convexmin(int n, double cumw[], double cs[], double y[])
{
    int	i, j, m;
    
    y[0] = cs[1]/cumw[1];
    for (i=2;i<=n;i++)
    {
        y[i-1] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
        if (y[i-2]>y[i-1])
        {
            j = i;
            while (y[j-2] > y[i-1] && j>1)
            {
                j--;
                if (j>1)
                    y[i-1] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
                else
                    y[i-1] = cs[i]/cumw[i];
                for (m=j;m<i;m++)	y[m-1] = y[i-1];
            }
        }
    }
}


double K(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y=(35.0/32)*pow(1-u,3);
    else
        y=0.0;
    
    return y;
}


double Kprime(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y = -(105.0/16)*x*pow(1-u,2);
    else
        y=0.0;
    
    return y;
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



void nelmin ( double fn ( double x[] ), int n, double start[], double xmin[],
             double *ynewlo, double reqmin, double step[], int konvge, int kcount,
             int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
    double ccoeff = 0.5;
    double del;
    double dn;
    double dnn;
    double ecoeff = 2.0;
    double eps = 0.001;
    int i;
    int ihi;
    int ilo;
    int j;
    int jcount;
    int l;
    int nn;
    double *p;
    double *p2star;
    double *pbar;
    double *pstar;
    double rcoeff = 1.0;
    double rq;
    double x;
    double *y;
    double y2star;
    double ylo;
    double ystar;
    double z;
    //
    //  Check the input parameters.
    //
    if ( reqmin <= 0.0 )
    {
        *ifault = 1;
        return;
    }
    
    if ( n < 1 )
    {
        *ifault = 1;
        return;
    }
    
    if ( konvge < 1 )
    {
        *ifault = 1;
        return;
    }
    
    p = new double[n*(n+1)];
    pstar = new double[n];
    p2star = new double[n];
    pbar = new double[n];
    y = new double[n+1];
    
    *icount = 0;
    *numres = 0;
    
    jcount = konvge;
    dn = ( double ) ( n );
    nn = n + 1;
    dnn = ( double ) ( nn );
    del = 1.0;
    rq = reqmin * dn;
    //
    //  Initial or restarted loop.
    //
    for ( ; ; )
    {
        for ( i = 0; i < n; i++ )
        {
            p[i+n*n] = start[i];
        }
        y[n] = fn ( start );
        *icount = *icount + 1;
        
        for ( j = 0; j < n; j++ )
        {
            x = start[j];
            start[j] = start[j] + step[j] * del;
            for ( i = 0; i < n; i++ )
            {
                p[i+j*n] = start[i];
            }
            y[j] = fn ( start );
            *icount = *icount + 1;
            start[j] = x;
        }
        //
        //  The simplex construction is complete.
        //
        //  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
        //  the vertex of the simplex to be replaced.
        //
        ylo = y[0];
        ilo = 0;
        
        for ( i = 1; i < nn; i++ )
        {
            if ( y[i] < ylo )
            {
                ylo = y[i];
                ilo = i;
            }
        }
        //
        //  Inner loop.
        //
        for ( ; ; )
        {
            if ( kcount <= *icount )
            {
                break;
            }
            *ynewlo = y[0];
            ihi = 0;
            
            for ( i = 1; i < nn; i++ )
            {
                if ( *ynewlo < y[i] )
                {
                    *ynewlo = y[i];
                    ihi = i;
                }
            }
            //
            //  Calculate PBAR, the centroid of the simplex vertices
            //  excepting the vertex with Y value YNEWLO.
            //
            for ( i = 0; i < n; i++ )
            {
                z = 0.0;
                for ( j = 0; j < nn; j++ )
                {
                    z = z + p[i+j*n];
                }
                z = z - p[i+ihi*n];
                pbar[i] = z / dn;
            }
            //
            //  Reflection through the centroid.
            //
            for ( i = 0; i < n; i++ )
            {
                pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
            }
            ystar = fn ( pstar );
            *icount = *icount + 1;
            //
            //  Successful reflection, so extension.
            //
            if ( ystar < ylo )
            {
                for ( i = 0; i < n; i++ )
                {
                    p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
                }
                y2star = fn ( p2star );
                *icount = *icount + 1;
                //
                //  Check extension.
                //
                if ( ystar < y2star )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = pstar[i];
                    }
                    y[ihi] = ystar;
                }
                //
                //  Retain extension or contraction.
                //
                else
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = p2star[i];
                    }
                    y[ihi] = y2star;
                }
            }
            //
            //  No extension.
            //
            else
            {
                l = 0;
                for ( i = 0; i < nn; i++ )
                {
                    if ( ystar < y[i] )
                    {
                        l = l + 1;
                    }
                }
                
                if ( 1 < l )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = pstar[i];
                    }
                    y[ihi] = ystar;
                }
                //
                //  Contraction on the Y(IHI) side of the centroid.
                //
                else if ( l == 0 )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
                    }
                    y2star = fn ( p2star );
                    *icount = *icount + 1;
                    //
                    //  Contract the whole simplex.
                    //
                    if ( y[ihi] < y2star )
                    {
                        for ( j = 0; j < nn; j++ )
                        {
                            for ( i = 0; i < n; i++ )
                            {
                                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                                xmin[i] = p[i+j*n];
                            }
                            y[j] = fn ( xmin );
                            *icount = *icount + 1;
                        }
                        ylo = y[0];
                        ilo = 0;
                        
                        for ( i = 1; i < nn; i++ )
                        {
                            if ( y[i] < ylo )
                            {
                                ylo = y[i];
                                ilo = i;
                            }
                        }
                        continue;
                    }
                    //
                    //  Retain contraction.
                    //
                    else
                    {
                        for ( i = 0; i < n; i++ )
                        {
                            p[i+ihi*n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                }
                //
                //  Contraction on the reflection side of the centroid.
                //
                else if ( l == 1 )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
                    }
                    y2star = fn ( p2star );
                    *icount = *icount + 1;
                    //
                    //  Retain reflection?
                    //
                    if ( y2star <= ystar )
                    {
                        for ( i = 0; i < n; i++ )
                        {
                            p[i+ihi*n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                    else
                    {
                        for ( i = 0; i < n; i++ )
                        {
                            p[i+ihi*n] = pstar[i];
                        }
                        y[ihi] = ystar;
                    }
                }
            }
            //
            //  Check if YLO improved.
            //
            if ( y[ihi] < ylo )
            {
                ylo = y[ihi];
                ilo = ihi;
            }
            jcount = jcount - 1;
            
            if ( 0 < jcount )
            {
                continue;
            }
            //
            //  Check to see if minimum reached.
            //
            if ( *icount <= kcount )
            {
                jcount = konvge;
                
                z = 0.0;
                for ( i = 0; i < nn; i++ )
                {
                    z = z + y[i];
                }
                x = z / dnn;
                
                z = 0.0;
                for ( i = 0; i < nn; i++ )
                {
                    z = z + pow ( y[i] - x, 2 );
                }
                
                if ( z <= rq )
                {
                    break;
                }
            }
        }
        //
        //  Factorial tests to check that YNEWLO is a local minimum.
        //
        for ( i = 0; i < n; i++ )
        {
            xmin[i] = p[i+ilo*n];
        }
        *ynewlo = y[ilo];
        
        if ( kcount < *icount )
        {
            *ifault = 2;
            break;
        }
        
        *ifault = 0;
        
        for ( i = 0; i < n; i++ )
        {
            del = step[i] * eps;
            xmin[i] = xmin[i] + del;
            z = fn ( xmin );
            *icount = *icount + 1;
            if ( z < *ynewlo )
            {
                *ifault = 2;
                break;
            }
            xmin[i] = xmin[i] - del - del;
            z = fn ( xmin );
            *icount = *icount + 1;
            if ( z < *ynewlo )
            {
                *ifault = 2;
                break;
            }
            xmin[i] = xmin[i] + del;
        }
        
        if ( *ifault == 0 )
        {
            break;
        }
        //
        //  Restart the procedure.
        //
        for ( i = 0; i < n; i++ )
        {
            start[i] = xmin[i];
        }
        del = eps;
        *numres = *numres + 1;
    }
    delete [] p;
    delete [] pstar;
    delete [] p2star;
    delete [] pbar;
    delete [] y;
    
    return;
}


void swap(double *x,double *y)
{
    double temp;
    temp=*x;
    *x=*y;
    *y=temp;
}


