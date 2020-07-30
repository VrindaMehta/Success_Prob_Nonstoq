#include "basic.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>
using namespace std;
typedef complex<double> dcomp;

#define lapack_complex_float complex<float>
#define lapack_complex_double complex<double>

#include <mkl_lapacke.h>


int main(int argc, char *argv[])    //argv[]=[readfile,t,dt,alpha,g]
{

    ::iota = -1;
    ::iota = sqrt(::iota);
    readfunc(argv[1]);
    int num_ev = 4;
    int i, j, step;
    double t, dt;
    vector<complex<double> > initialgroundstate(p), finalgroundstate(p);
    t = stod(argv[2]);
    dt = stod(argv[3]);
    cout<<endl;
    step = t/dt;
    cout<<"#t =  "<<t<<endl;
    cout<<"#dt = "<<dt<<endl;
    cout<<"#Number of time steps = "<<step<<endl;
    cout<<"# g = "<<argv[5]<<endl<<endl;
    cout<<"#Trigger = "<<argv[4]<<endl;
    cout<<endl;

    results_spectrum(argv[4], dt, num_ev, step, stof(argv[5]), finalgroundstate);
    results_productevolution(dt, step, finalgroundstate);

    delete [] h_x;
    delete [] h1_x;
    delete [] h_y;
    delete [] h_z;
    delete [] J_x;
    delete [] J_y;
    delete [] J_z;
    delete [] flag;
    return 0;

}
