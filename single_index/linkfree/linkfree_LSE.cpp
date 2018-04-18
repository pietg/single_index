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
#include <random>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

int n;
double **xx,**S,*yy,*vv,*f,*psi;

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
void    gen_data(int m, int n, double alpha[], double **xx, double yy[], double vv[], int seed);
double  best_nearby (int m, double delta[], double point[], double prevbest,
                     double f(int m, double alpha[]), int *funevals);
int     hooke(int m, double startpt[], double endpt[], double rho,
              double eps, int itermax, double f(int m, double alpha[]));
int     CompareTime(const void *a, const void *b);
void    statistics(int n, int m, double **estimate, double *mean, double **covar);


// [[Rcpp::export]]

List ComputeLinkFree_LSE()
{
    int             i,j,m,iter,iter2,NumIt,seed, nIterations;
    double          *alpha,*alpha0,f2, rho;
    double          **covar,*mean,**estimate, *alpha_init;
    double          eps,sum;
    clock_t         StartTime, StopTime;
    double          Time_simulation;
    
    Rcout << std::endl << std::endl;
    Rcout << "Piet Groeneboom 2018" << std::endl << "For more information see:" << std::endl;
    Rcout << "Score estimation in the monotone single index model" << std::endl;
    Rcout << "Fadoua Balabdaoui, Piet Groeneboom and Kim Hendrickx" << std::endl;
    Rcout << "https://arxiv.org/abs/1712.05593" << std::endl << std::endl;
    
    seed=1152;
    
    n=1000;
    NumIt=1000;
    m=3;
    
    Rcout << "Number of observations:" << std::setw(10) << n << std::endl << std::endl;
    Rcout << "Number of samples:" << std::setw(10) << NumIt << std::endl << std::endl;
    
    nIterations=100;
    
    iter2=1000;
    eps=1.0e-10;
    
    xx = new double *[n];
    for (i=0;i<n;i++)
        xx[i] = new double [m];
    
    S = new double *[m];
    for (i=0;i<m;i++)
        S[i] = new double [m];
    
    estimate = new double *[NumIt];
    for (i=0;i<NumIt;i++)
        estimate[i] = new double [m];
    
    covar = new double *[m];
    for (i=0;i<m;i++)
        covar[i] = new double [m];
    
    mean = new double[m];
    
    yy= new double[n];
    vv= new double[n];
    alpha= new double[m];
    alpha0= new double[m];
    alpha_init= new double[m];
    f= new double[m];
    psi  = new double[n];
    
    for (i=0;i<m;i++)
        alpha0[i]=alpha[i]=alpha_init[i]=1/sqrt(m);
    
    
    StartTime = clock();

    
    Rcout << "     Iteration  " << "  alpha1  " << "      alpha2  "<< "       alpha3  " << "        criterion  "
    << "             norm  "  <<  std::endl << std::endl;

    
    for (iter=0;iter<NumIt;iter++)
    {
        gen_data(m,n,alpha0,xx,yy,vv,seed+iter);
        
        for (i=0;i<m;i++)
            alpha_init[i]=1/sqrt(m);
    
        rho=0.5;
        
        iter2=hooke(m,alpha_init,alpha,rho,eps,nIterations,&criterion);
        
        for (j=0;j<m;j++)
            estimate[iter][j]=alpha[j];
        
        f2 = criterion(m,alpha);
        
        f2=sqrt(f2);
        
        sum=0;
        
        for (i=0;i<m;i++)
        {
            for (j=0;j<m;j++)
                sum += alpha[i]*S[i][j]*alpha[j];
        }
        
        sum=sqrt(sum);
        
        //printf("%5d   %15.10f   %15.10f   %15.10f   %15.10f     %15.10f\n",iter+1,alpha[0],alpha[1],alpha[2],f2,sum);
                
        Rcout  << setw(10) << iter+1 << setprecision(6) <<  setw(15) << alpha[0] << setprecision(6) <<  setw(15) << alpha[1] << setprecision(6) <<  setw(15) << alpha[2] <<  setprecision(10) << setw(20) << f2 << setprecision(10) << setw(15) << sum << std::endl;
    }
    
    StopTime  = clock();
    Time_simulation   = (double)(StopTime - StartTime)/(double)CLOCKS_PER_SEC;
    
    Rcout << std::endl << std::endl;
    Rcout << "The computations took    " << setprecision(10) << Time_simulation << "   seconds"  << std::endl;
    
    statistics(NumIt,m,estimate,mean,covar);

    NumericMatrix out1 = NumericMatrix(NumIt,m);
    
    for (i=0;i<NumIt;i++)
    {
        for (j=0;j<m;j++)
            out1(i,j) = estimate[i][j];
    }

    NumericVector out2 = NumericVector(m);
    
    for (i=0;i<m;i++)
        out2(i)=mean[i];
    
    
    NumericMatrix out3 = NumericMatrix(m,m);
    
    for (i=0;i<m;i++)
    {
        for (j=0;j<m;j++)
            out3(i,j)=n*covar[i][j];
    }
    
    
    ofstream file0_("estimate.txt");
    
    if (file0_.is_open())
    {
        for (i=0;i<NumIt;i++)
        {
            for (j=0;j<m;j++)
                file0_ << setprecision(11) << setw(20) << estimate[i][j];
            file0_ << "\n";
        }
        file0_.close();
    }
    

    Rcout << "Making output list" << std::endl;
    
    // make the list for the output
    
    List out = List::create(Rcpp::Named("LFLSE")=out1,Rcpp::Named("means")=out2,Rcpp::Named("covariance_matrix")=out3);
    
    
    Rcout << "Freeing memory" << std::endl;
    
    // free memory
    
    delete[] mean, delete[] yy, delete[] vv, delete[] alpha0, delete[] alpha,
    delete[] alpha_init, delete[] f, delete[] psi;
    
    for (i = 0;i < n;i++) delete[] xx[i];
    delete[] xx;
    
    for (i = 0;i < m;i++) delete[] S[i];
    delete[] S;
    
    for (i = 0;i < NumIt;i++) delete[] estimate[i];
    delete[] estimate;
    
    for (i = 0;i < m;i++) delete[] covar[i];
    delete[] covar;
    
    return out;
    
}

void gen_data(int m, int n, double alpha[], double **xx, double yy[], double vv[], int seed)
{
    int	i,j;
    
    
    std::mt19937 generator(seed);
    std::normal_distribution<double> dis_normal(0.0,1.0);
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<m;j++)
            xx[i][j] = dis_normal(generator);
        
        yy[i]=0;
        
        for (j=0;j<m;j++)
            yy[i] += alpha[j]*xx[i][j];
        
        yy[i] = pow(yy[i],3)+dis_normal(generator);
    }
    
    sort_alpha(m,n,xx,alpha,vv,yy);
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

double criterion(int m, double alpha[])
{
    int i,j,k;
    double lambda,sum,*x_mean,*g;
    
    x_mean  = new double[m];
    g       = new double[m];
    
    for (j=0;j<m;j++)
    {
        sum=0;
        for (i=0;i<n;i++)
            sum += xx[i][j];
        x_mean[j] = sum/n;
    }
    
    
    for (i=0;i<m;i++)
    {
        for (j=0;j<m;j++)
            S[i][j]=0;
    }
    
    for (j=0;j<m;j++)
    {
        for (k=0;k<m;k++)
        {
            for (i=0;i<n;i++)
                S[j][k]+= (xx[i][j]-x_mean[j])*(xx[i][k]-x_mean[k]);
            S[j][k]/=n;
        }
    }
    
    for (j=0;j<n;j++)
        psi[j]=0;
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<m;j++)
            psi[i]+=alpha[j]*(xx[i][j]-x_mean[j]);
    }
    
    
    for (j=0;j<m;j++)
        f[j]=g[j]=0;
    
    for (j=0;j<m;j++)
    {
        for (i=0;i<n;i++)
            f[j] += (xx[i][j]-x_mean[j])*(yy[i]-psi[i]);
        f[j]/=n;
    }
    
    lambda=0;
    
    for (i=0;i<m;i++)
        lambda += alpha[i]*f[i];
    
    for (i=0;i<m;i++)
    {
        for (k=0;k<m;k++)
            g[i] += S[i][k]*alpha[k];
    }
    
    for (i=0;i<m;i++)
        f[i] -= lambda*g[i];
    
    sum=0;
    
    for (i=0;i<m;i++)
        sum += SQR(f[i]);
    
    delete[] x_mean;
    delete[] g;
    
    return sum;
}


void statistics(int n, int m, double **estimate, double *mean, double **covar)
{
    int i,j,k;
    double sum;
    
    for (j=0;j<m;j++)
    {
        sum=0;
        for (i=0;i<n;i++)
            sum += estimate[i][j];
        mean[j]=sum/n;
    }
    
    for (j=0;j<m;j++)
    {
        for (k=0;k<m;k++)
        {
            sum=0;
            for (i=0;i<n;i++)
                sum += (estimate[i][j]-mean[j])*(estimate[i][k]-mean[k]);
            covar[j][k] = sum/n;
        }
    }
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



int CompareTime(const void *a, const void *b)
{
    if ((*(data_object *) a).v < (*(data_object *) b).v)
        return -1;
    if ((*(data_object *) a).v > (*(data_object *) b).v)
        return 1;
    return 0;
}


