#include "basic.h"
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "omp.h"
#include <random>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <vector>

using namespace std;

#define lapack_complex_float complex<float>
#define lapack_complex_double complex<double>

#include <mkl_lapacke.h>

void Print(vector<complex<double> > &psi) { // Prints the wavefunction

    cout<<"The wavefunction is "<<endl;
    for(int i = 0; i < p; i++) {
        cout<<"#"<<psi[i]<<endl;
    }
}

double Norm(vector<complex<double> > &psi, char r) { //Returns the norm of the wavefunction. Normalises it if response r='y'

    double norm = 0;

    for(int i = 0; i < p; i++){
        norm += pow( real(psi[i]),2 ) + pow( imag(psi[i]),2 );
    }
    norm = sqrt(norm);
    if(norm>pow(10,-10)) {
        if (r == 'y') {
            for(int i = 0; i < p; i++){
                psi[i] /= norm;
            }
        }
    }
    return norm;
}


int X_(vector<complex<double> >& psi) {

    int i, l1, l2, k1, k0, a;
	for (i = 0; i < n; i++) {
        a = 1 << i;

    for (k1 = 0; k1 < p; k1 += (2 * a)){
        for (k0 = 0; k0 < a; k0++) {
                l1 = (k0 | k1);
                l2 = (l1 | a);
                vector<complex<double> > vec = {psi[l1],psi[l2]};
                psi[l1] = (vec[0] - ::iota * vec[1]) / sqrt(2);
                psi[l2] = (-::iota* vec[0] + vec[1]) / sqrt(2);

            }
        }

    }
	return 0;
}

int X(vector<complex<double> >& psi) {

    int i, l1, l2, k1, k0, a;
	for (i = 0; i < n; i++){                //X
        a = 1 << i;
        for (k1 = 0; k1 < p; k1 += (2 * a)) {
           for (k0 = 0; k0 < a; k0++) {
                l1 = (k0 | k1);
                l2 = (l1 |a);
                vector<complex<double> > vec = {psi[l1],psi[l2]};
                psi[l1] = (vec[0] + ::iota * vec[1]) / sqrt(2);
                psi[l2] = (::iota * vec[0] + vec[1]) / sqrt(2);

            }
        }
    }

	return 0;
}


int Y(vector<complex<double> >& psi) {

    int i, l1, l2, k1, k0, a;
	for ( i = 0; i < n; i++){
        a = 1 << i;
        for (k1 = 0; k1 < p; k1 += (2 * a)) {
            for (k0 = 0; k0 < a; k0++) {
                l1 = (k0 | k1);
                l2 = (l1 | a);
                vector<complex<double> > vec = {psi[l1],psi[l2]};
                psi[l1] = (vec[0] + vec[1]) / sqrt(2);
                psi[l2] = (-vec[0] + vec[1]) / sqrt(2);
            }
        }
    }

    return 0;
}



int Y_(vector<complex<double> >& psi) {

    int i, l1, l2, k1, k0, a;
	for (i = 0; i < n; i++) {
        a = 1 << i;
        for (k1 = 0; k1 < p; k1 += (2 * a)) {
            for (k0 = 0; k0 < a; k0++) {
                l1 = (k0 | k1);
                l2 = (l1 | a);
                vector<complex<double> > vec = {psi[l1],psi[l2]};
                psi[l1] = (vec[0] - vec[1]) / sqrt(2);
                psi[l2] = (vec[0] + vec[1]) / sqrt(2);
            }
        }
    }

	return 0;
}

void successprob(vector<complex<double> > psi3, vector<complex<double> > groundstate) {
  complex<double> prob = 0;
  for (int i = 0; i < p; i++) {
    prob +=conj(psi3[i]) * groundstate[i];
  }
  prob = pow(real(prob), 2) + pow(imag(prob), 2);
  cout<<setw(20)<<real(prob);
}

void Magnetization(vector<complex<double> > psi3) {
  double sum;
  int a;
  for (int i = 0; i < n; i++) {
    sum = 0;
    a = 1 << i;
    for (int k = 0; k < p; k++) {
      if( (a & k) > 0) { //checking if ith bit is not zero
        sum += pow(real(psi3[k]),2) + pow(imag(psi3[k]),2);
      }
      else {
        sum -= pow(real(psi3[k]),2) + pow(imag(psi3[k]),2);
      }
    }
    cout<<setw(20)<<-sum;
  }
}

void Energy(vector<complex<double> > psi, double s) {

    vector<complex<double> > psi2(p);
    double energy = 0;

    Hamil(psi,psi2,s);

    for (int i=0;i<p;i++) {
        energy+=real(conj(psi[i])*psi2[i]);
    }
    cout<<setw(20)<<energy;

}
