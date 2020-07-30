#include <vector>
#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include "basic.h"
using namespace std;

#define lapack_complex_float complex<float>
#define lapack_complex_double complex<double>

#include <mkl_lapacke.h>


void Product_Single(vector<complex<double> >& psi, double dt,double s)
{
    int i,l1,l2,k1,k0,a;
	dcomp Mul[2][2];
    double h_i;
    for (i=0;i<n;i++)
    {
        h_i= pow((pow(((1-s)*h_x[i]+s*(1-s)*h1_x[i]),2)+ pow(((1-s)*s*h_y[i]),2)+ pow((s*h_z[i]),2)),0.5);
        a=1<<i;
        if (h_i!=0)
        {
            Mul[0][0]= cos(h_i*dt)+ (::iota*s*h_z[i]*sin(h_i*dt)/h_i);
            Mul[0][1]= (::iota*((1-s)*h_x[i]+s*(1-s)*h1_x[i])+ s*(1-s)*h_y[i])*sin(h_i*dt)/h_i;
            Mul[1][0]= (::iota*((1-s)*h_x[i]+s*(1-s)*h1_x[i])- s*(1-s)*h_y[i])*sin(h_i*dt)/h_i;
            Mul[1][1]= cos(h_i*dt)- (::iota*s*h_z[i]*sin(h_i*dt)/h_i);
            for (k1=0;k1<p;k1+=(2*a))
            {
                for (k0=0;k0<a;k0++)
                {
                    l1=(k0|k1);
                    l2=(l1|a);
                    dcomp vec[2]={psi[l1],psi[l2]};
                    psi[l1]= Mul[0][0]*vec[0]+ Mul[0][1]*vec[1];
                    psi[l2]= Mul[1][0]*vec[0]+ Mul[1][1]*vec[1];
                }
            }
        }
    }
}

void Product_Z(vector<complex<double> >& psi1,double dt,double s)
{
    int i,j,k,a,b;
    double sum;
    for (k=0;k<p;k++)
    {
        sum=0;
        for (i=0;i<n;i++)
        {
            a=1<<i;
            for (j=i+1;j<n;j++)
            {
                b=1<<j;
                if (((a&k)>>i)==((b&k)>>j))
                    sum-=s*J_z[j+i*n];
                else
                    sum+=s*J_z[j+i*n];
            }
        }
        psi1[k]*=exp(-::iota*dt*sum);
    }
}


void Product_H_X(vector<complex<double> >& psi, double dt,double s)
{
    int i,j,k,a,b;
    double sum;

    Y(psi);

    for (k=0;k<p;k++)
    {
    sum=0;
    for (i=0;i<n;i++)
    {


        a=1<<i;
        for (j=i+1;j<n;j++)
        {
            b=1<<j;
            if (((a&k)>>i)==((b&k)>>j))
                    sum-=s*(1-s)*J_x[j+i*n];
                else
                    sum+=s*(1-s)*J_x[j+i*n];
        }

    }
        psi[k]*=exp(-::iota*dt*sum);
}


    Y_(psi);
}


void Product_H_Y(vector<complex<double> >& psi, double dt,double s)
{
    int i,j,k,a,b;
    double sum;


    X_(psi);
    for (k=0;k<p;k++)
    {
    sum=0;
    for (i=0;i<n;i++)
    {
        a=1<<i;
        for (j=i+1;j<n;j++)
        {
            b=1<<j;
            if (((a&k)>>i)==((b&k)>>j))
                    sum-=s*J_y[j+i*n];
                else
                    sum+=s*J_y[j+i*n];
        }

    }
        psi[k]*=exp(-::iota*dt*sum);
    }

    X(psi);

}


void productevolution(vector<complex<double> >& psi, double dt,double s)
{
    Product_Z(psi,dt/2,s);
    Product_H_Y(psi,dt/2,s);
    Product_H_X(psi,dt/2,s);
    Product_Single(psi,dt,s);
    Product_H_X(psi,dt/2,s);
    Product_H_Y(psi,dt/2,s);
    Product_Z(psi,dt/2,s);
}
