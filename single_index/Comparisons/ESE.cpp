//
//  link-free estimate in the single index model
//
//  Created by Piet Groeneboom on 16/04/18.
//  Copyright (c) 2018 Piet Groeneboom. All rights reserved.
//
//  The Hooke-Jeeves algorithm of the TOMS178 library is used:
//    The ALGOL original is by Arthur Kaupe.
//    C version by Mark Johnson
//    C++ version by John Burkardt

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
//#include <chrono>
#include <random>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

int n;
double **xx,*yy,*vv,*cumw,*cs,*f,*psi,*derivative;

#define SQR(x) ((x)*(x))

typedef struct
{
    int index;
    double v;
    double y;
}
data_object;


void    swap(double *x,double *y);
double  criterion(int m, double alpha[]);
void    sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[]);
void    convexmin(int n, double cumw[], double cs[], double y[]);
double  best_nearby (int m, double delta[], double point[], double prevbest,
                     double f(int m, double alpha[]), int *funevals);
int     hooke(int m, double startpt[], double endpt[], double rho, double eps,
              int itermax, double f(int m, double alpha[]));
int     CompareTime(const void *a, const void *b);
double  KK(double x);
double  K(double x);
double  Kprime(double x);


// [[Rcpp::export]]

List ComputeESE(NumericMatrix X, NumericVector y, NumericVector alpha0, int m1)
{
    int             i,j,m,iter,nIterations;
    double          *alpha,rho,f2;
    double          *alpha_init;
    double          eps;
    
    // determine the sample size
    
    n = (int)(y.size());
    
    // m is the dimension
    
    m= (int)m1;
    
    nIterations=1000;

    eps=1.0e-5;
    rho=0.5;
    
    xx = new double *[n];
    for (i=0;i<n;i++)
        xx[i] = new double [m];
    
    yy= new double[n];
    vv= new double[n];
    alpha= new double[m];
    alpha_init= new double[m];
    f= new double[m];
    psi  = new double[n+1];
    derivative  = new double[n+1];
    
    cumw = new double[n+1];
    cs = new double[n+1];
    
    cumw[0]=cs[0]=0;
    
    for (i=0;i<n;i++)
        yy[i]=(double)y[i];
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<m;j++)
            xx[i][j]=(double)X(i,j);
    }
    
    for (i=0;i<m;i++)
        alpha_init[i]=(double)alpha0(i);
    
    iter=hooke(m,alpha_init,alpha,rho,eps,nIterations,&criterion);
    
    f2 = criterion(m,alpha);
    
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

    // Number of iterations of Hooke-Jeeves
    
    int out4 = iter;
    
    // make the list for the output, containing alpha and the estimate of psi
    
    List out = List::create(Rcpp::Named("data")=out0,Rcpp::Named("alpha")=out1,Rcpp::Named("psi")=out2,
                                        Rcpp::Named("derivative")=out3,Rcpp::Named("niter")=out4);
    
    // The last element shows the number of iterations in the Hooke-Jeeves algorithm
    
    // free memory
    
    delete[] yy, delete[] vv, delete[] alpha,
    delete[] alpha_init, delete[] f, delete[] psi, delete[] cumw, delete[] cs,
    delete[] derivative;
    
    for (i = 0;i < n;i++) delete[] xx[i];
    delete[] xx;
    
    return out;
    
}

double criterion(int m, double alpha[])
{
    int i,j,m1;
    double h1,h,lambda,sum,*uu,*pp,A,B;
    
    uu= new double[n];
    pp= new double[n];
    
    sum=0;
    
    for (i=0;i<m;i++)
        sum += SQR(alpha[i]);
    
    for (i=0;i<m;i++)
        alpha[i]/=sqrt(sum);
    
    sort_alpha(m,n,xx,alpha,vv,yy);
    
    A=vv[0];
    B=vv[n-1];
    
    h1=0.5*(B-A)*pow((double)n,-1.0/5);
    h=0.5*(B-A)*pow((double)n,-1.0/7);
    
    //printf("%15.10f\n",h);
    
    for (i=1;i<=n;i++)
    {
        cumw[i]=i*1.0;
        cs[i]=cs[i-1]+yy[i-1];
    }
    
    convexmin(n,cumw,cs,psi);
    
    j=-1;
    
    for (i=1;i<n;i++)
    {
        if (psi[i]>psi[i-1])
        {
            j++;
            uu[j]=vv[i];
            pp[j]=psi[i]-psi[i-1];
        }
    }
    
    m1=j;
    
    for (i=0;i<n;i++)
    {
        if (vv[i]<=B-h && vv[i]>=A+h)
        {
            sum=0;
            for (j=0;j<m1;j++)
                sum+= K((vv[i]-uu[j])/h)*pp[j]/h;
            derivative[i]=sum;
        }
        
        if (vv[i]>B-h)
        {
            sum=0;
            for (j=0;j<m1;j++)
            {
                sum += K((B-h-uu[j])/h)*pp[j]/h;
                sum += (vv[i]-B+h)*Kprime((B-h-uu[j])/h)*pp[j]/SQR(h);
            }
            derivative[i]=sum;
        }
        
        if (vv[i]<A+h)
        {
            sum=0;
            for (j=0;j<m1;j++)
            {
                sum += K((A+h-uu[j])/h)*pp[j]/h;
                sum += (vv[i]-A-h)*Kprime((A+h-uu[j])/h)*pp[j]/SQR(h);
            }
            derivative[i]=sum;
        }
    }
    
    
    for (j=0;j<m;j++)
        f[j]=0;
    
    for (j=0;j<m;j++)
    {
        for (i=0;i<n;i++)
            f[j] += xx[i][j]*derivative[i]*(psi[i]-yy[i]);
    }
    
    for (j=0;j<m;j++)
        f[j]/=n;
    
    lambda=0;
    
    for (i=0;i<m;i++)
        lambda += alpha[i]*f[i];
    
    for (i=0;i<m;i++)
        f[i] -= lambda*alpha[i];
    
    sum=0;
    
    for (i=0;i<m;i++)
        sum += SQR(f[i]);
    
    delete[] uu; delete[] pp;
    
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



double KK(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y = (16.0 + 35*x - 35*pow(x,3) + 21*pow(x,5) - 5*pow(x,7))/32.0;
    else
    {
        if (x>1)
            y=1;
        else
            y=0;
        
    }
    
    return y;
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
    
    // printf("\n %16.10f %15.10f %15.10f\n \n", newx[0],newx[1],newx[2]);
    
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


