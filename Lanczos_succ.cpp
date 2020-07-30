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


void Initialise(vector<complex<double> > &psi, int seed) { // Randomly initialise the wavefunction based on the seed. Every coefficient [-0.5,0.5]

    srand(seed);

    for(int i = 0; i < p;i++){
        psi[i] = { 2 * ( double(rand()) / RAND_MAX) - 1, 2 * ( double(rand()) / RAND_MAX) - 1 };
    }
}

void Hamil_Z(vector<complex<double> >& psi1,double s) {   //Alters psi to Hz(psi). Hz = sum_i(h_z^i. sigma_z^i) - sum_i,j(J_z^i,j. sigma_z^i. sigma_z^j)

    int i, j, k, a, b;
    double sum;
    for (k = 0; k < p; k++) {
        sum = 0;
        for (i = 0; i < n; i++) {
            a = 1 << i;
            if( (a & k) > 0) { //checking if ith bit is not zero
                sum += s * h_z[i];
            }
            else{
                sum -= s * h_z[i];
            }
            for (j = i+1; j < n; j++){
                b = 1 << j;
                if ( ((a & k) >> i) == ( (b & k) >> j)) {
                    sum -= s * J_z[j + i *n];
                }
                else {
                    sum += s * J_z[j + i * n];
                }
            }
        }
        psi1[k] *= sum;
    }
}



void Hamil_X(vector<complex<double> >& psi, double s) {

    int i, j, k, a, b;
    double sum;

    Y(psi);

    for (k = 0; k < p; k++) {
        sum = 0;

        for (i = 0; i < n; i++) {
            a = 1 << i;
            if( (a & k) > 0){ //checking if ith bit is not zero
                sum += (1 - s) * h_x[i] + s * (1 - s) * h1_x[i];
            }
            else{
                sum -= (1 - s) * h_x[i] + s * (1 - s) * h1_x[i];
            }
            for (j = i+1; j < n; j++){
                b = 1 << j;
                if ( ( (a & k) >> i) == ( (b & k) >> j) ){
                    sum -= s * (1 - s) * J_x[j + i * n];
                }
                else{
                    sum += s * (1 - s) * J_x[j + i * n];
                }
            }
        }

        psi[k] *= sum;
    }

    Y_(psi);
}


void Hamil_Y(vector<complex<double> >& psi, double s) {

    int i, j, k, a, b;
    double sum;

    X_(psi);
    for (k = 0; k < p; k++) {
        sum = 0;
        for (i = 0; i < n; i++) {
            a = 1 << i;
            if ( (a & k) > 0) { //checking if ith bit is not zero
                sum += s * h_y[i];
            }
            else{
                sum -= s * h_y[i];
            }
            for (j = i + 1; j < n; j++){
                b = 1 << j;
                if ( ( (a & k) >> i) == ( (b & k) >> j) ) {
                 sum -= s * J_y[j + i * n];
                }
                else{
                sum += s * J_y[j + i * n];
                }
            }
        }
        psi[k] *= sum;
    }

    X(psi);
}

void Hamil(vector<complex<double> > psi, vector<complex<double> > &psi2, double s){ //Alters psi2 to psi2 = (Hx + Hy + Hz) psi

    vector<complex<double> > psi1(p); //pseudo vector for storing the original wave vector
    int i;

    psi1 = psi;

    Hamil_Z(psi1, s);

    psi2 = psi1;
    psi1 = psi;

    Hamil_X(psi1, s);

    for(i = 0; i < p; i++){
        psi2[i] += psi1[i];
    }
    psi1 = psi;

    Hamil_Y(psi1, s);

    for (i = 0; i < p; i++){
        psi2[i] += psi1[i];  //psi2 contains the final result. psi unchanged.
    }
}


void inst_overlap(int num_ev, vector<complex<double> > psi3, vector<complex<double> > Eigenvect ) {

    complex<double> o0 = 0, o1 = 0, o2 = 0;

    for (int i = 0; i < p; i++) {
        o0 += (psi3[i]) * conj( Eigenvect[i * num_ev] );
        o1 += conj(psi3[i]) * Eigenvect[1 + i * num_ev];
        o2 += conj(psi3[i]) * Eigenvect[2 + i * num_ev];
    }

    o0 = pow (real(o0),2) + pow (imag(o0),2);
    o1 = pow (real(o1),2) + pow (imag(o1),2);
    o2 = pow (real(o2),2) + pow (imag(o2),2);

    cout<<setw(20)<<real(o0)<<setw(20)<<real(o1)<<setw(20)<<real(o2);
}

void StandardDeviation(vector<complex<double> > Eigenvect, double s, int num_ev) {

    double E_sq, energy, std;
    vector<complex<double> > psi1(p), psi2(p), psi3(p);

    for (int i = 0; i < num_ev; i++) {
        for (int j = 0; j < p; j++)
            psi1[j] = Eigenvect[i + j * num_ev];
        energy = 0;

        Hamil(psi1, psi2, s);

        for (int j = 0; j < p; j++) {
            energy += real( conj(psi1[j]) * psi2[j]);
        }

        for (int j = 0; j < p; j++){
            psi3[j] = psi2[j] - energy * psi1[j];
        }

        std = Norm(psi3,'n');
        cout<<setw(20)<<std;
    }
}

void Egev(int k, vector<double> diag, vector<double> subdiag, vector<complex<double> > V, vector<double> &eigenvalue){

    int i, j;
    int info = LAPACKE_zsteqr(LAPACK_ROW_MAJOR,'N', k, &diag[0], &subdiag[0], &V[0], k); //give 1 as another argument (for LDZ) if this gives error.
    if (info != 0)
        cout<<"Eigenvalues cannot be computed"<<endl;
    for(i = 0; i < k; i++)
        eigenvalue[i]=diag[i];
}

void Egvector(int obtained_k, int num_ev, vector<double> diag, vector<double> subdiag, vector<double> eigenvalue, vector<complex<double> > &V){
    vector<int> iblock(obtained_k,1), isplit(obtained_k,(obtained_k+1)), ifail(num_ev); //num_ev is the desired number, but k eigenvalues
    isplit[0] = obtained_k;

    int info = LAPACKE_zstein(LAPACK_ROW_MAJOR, obtained_k, &diag[0], &subdiag[0], num_ev, &eigenvalue[0], &iblock[0], &isplit[0], &V[0],num_ev, &ifail[0]);  //remove Z if this gives an error. Remove also the declaration in Lanczos().
     if (info != 0)
        cout<<"Eigenvectors cannot be computed. Info = "<<setw(10)<<info<<endl;
}



void Lanczos(vector<complex<double> >& Eigenvect, double s, int num_ev) {
    int i, j, k, l, obtained_k;
    int k_max = 160, lookback = 160;
    double diff , a = 0.0, b = 0.0;
    complex<double>  overlap;
    vector<double>  eigenvalue(k_max), old_eigenvalue(k_max, 0.0), diag(k_max), subdiag(k_max-1);
    vector<complex<double> > psi(p), psi1(p), psi2(p,0), U(p * k_max), V(k_max * num_ev), test_vect(num_ev * p, 0); //psi1 pseudo vector so that it can be normalised.

    Initialise(psi , 4);
    Norm(psi, 'y');

    // Fill the eigenvect with zeroes
    fill(Eigenvect.begin(), Eigenvect.end(), 0);

    for(int j = 0; j < p; j++){
        U[j * k_max] = psi[j];
    }

    psi1 = psi;             //so that psi remains unchanged and can be subsequently used for product evolution

    Hamil(psi1, psi2, s);

    //#pragma omp parallel for simd reduction(+:a)
    for(j = 0; j < p; j++){
        a += real(conj(psi1[j]) * psi2[j]);
    }

    diag[0] = a;

    for(k = 1; k < k_max; k++)     //if k_max less than num_ev==> error
    {
       //cout<<k<<setw(20)<<a<<setw(20)<<b<<endl;
        if(k == 1) {
            //#pragma omp parallel for
            for(j = 0; j < p; j++) {
                psi1[j] = psi2[j] - a * psi1[j];
            }
        }

        else {
            //#pragma omp parallel for
            for(j = 0; j < p; j++){
                psi1[j] = psi2[j] - a * psi1[j] - b * U[(k-2) + j * k_max];
            }
        }

        for(l = 0; l < k; l++) {
            overlap = 0;
            for(j = 0; j < p; j++) {
                  overlap+=conj(U[l+j*k_max])*psi1[j];
            }
            for(j = 0; j < p; j++){
                psi1[j] -= overlap * U[l + j * k_max];
            }
        }

        b=Norm(psi1,'y');

        if(b < pow(10,-10) ) {
          Initialise(psi1, k);

          for(l = 0; l < k; l++) {
              overlap = 0;
              for(j = 0; j < p; j++) {
                  overlap += conj(U[l + j * k_max]) * psi1[j];
              }
              for(j = 0; j < p; j++){
                  psi1[j] -= overlap * U[l + j * k_max];
              }

              Norm(psi1,'y');
          }
          //b=Norm(psi1,'y');
        }

        for(j=0;j<p;j++){
            U[k+j*k_max]=psi1[j];
        }

        Hamil(psi1, psi2, s);

        a = 0.0;
        // #pragma omp parallel for simd reduction(+:a)
        for(j = 0; j < p; j++) {
            a += real( conj(psi1[j]) * psi2[j] );     //a_i
        }

        diag[k] = a;
        subdiag[k-1] = b;

        if(k >= num_ev)
        {
            Egev(k+1, diag, subdiag, V, eigenvalue);
            diff = 0;
            for(j = 0; j < num_ev; j++) {
                diff += abs(old_eigenvalue[j] - eigenvalue[j]);
            }
                if(k % lookback == 0 && diff<pow(10,-10)) {
                    break;
                }

                else {
                    for(j = 0; j < num_ev; j++) {
                        old_eigenvalue[j] = eigenvalue[j];
                    }
                }
        }
    }
    obtained_k = k;

    Egvector(obtained_k, num_ev, diag, subdiag, eigenvalue, V);

    for(l = 0; l < p; l++) {
        for(i = 0; i < num_ev; i++) {
            for(j = 0; j < obtained_k; j++) {
		                  Eigenvect.at(i + l * num_ev) += U[j + l * k_max] * V[i + j * num_ev];
		        }
        }
    }


    for(int j = 0; j < num_ev; j++)
        cout<<setw(20)<<eigenvalue[j];

    StandardDeviation(Eigenvect, s, num_ev);
    cout<<setw(20)<<obtained_k;


}

void Lanczos_evolution(vector<complex<double> >& U, vector<complex<double> >& V, vector<double>& eigenvalue, vector<complex<double> >& psi, double s, int obtained_k, int num_ev, double dt)
{

    int i,j;
    vector<complex<double> > psi1(p,0);  //pseudo
    for(i=0;i<obtained_k;i++)
    {
        for(j=0;j<p;j++)
            psi1[i]+=conj(U[j+i*p])*psi[j];
    }

    for(i=0;i<obtained_k;i++)
    {
        psi[i]=psi1[i];
        psi1[i]=0;
    }
    for(i=0;i<num_ev;i++)
    {
        for(j=0;j<obtained_k;j++)
            psi1[i]+=conj(V[i+j*num_ev])*psi[j];
    }
    for(i=0;i<obtained_k;i++)
    {
        psi[i]=psi1[i];
        psi1[i]=0;

    }
    for(i=0;i<num_ev;i++)
        psi[i]*=exp(-::iota*eigenvalue[i]*dt);
    for(i=0;i<obtained_k;i++)
    {
        for(j=0;j<num_ev;j++)
            psi1[i]+=V[j+i*num_ev]*psi[j];
    }
    for(i=0;i<obtained_k;i++)
    {
        psi[i]=psi1[i];
        psi1[i]=0;
    }
    for(i=0;i<p;i++)
    {
        for(j=0;j<obtained_k;j++)
            psi1[i]+=U[i+j*p];
    }

    for(i=0;i<p;i++)
        psi[i]=psi1[i];

}

void results_spectrum(string alpha, double dt, int num_ev, int step, float g, vector<complex<double> >& finalgroundstate)
{

    int i;
    double s;
    vector<complex<double> > Eigenvect(p * num_ev), psi3(p,1);
    Norm(psi3,'y');

    cout<<setprecision(4);
    cout<<"# s p instantaneous_Energy overlap_1,2,3 E_1....10 StandardDeviation_1...10"<<endl<<endl;
    cout<<std::fixed;
    
    params(alpha,g);
    
    s = 1;
    
    cout<<setprecision(4)<<scientific<<s;
    cout<<setprecision(10);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    Lanczos(Eigenvect,s,num_ev);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>( end - start ).count();
    cout<<setw(20)<<"time:"<<diff/1e3;
    
    for (int j = 0; j < p; j++) {
        finalgroundstate[j] = Eigenvect[j * num_ev];
    }
    
    cout<<endl;

}

void results_productevolution(double dt, int step, vector<complex<double> > finalgroundstate ) {
  double s, norm;
  vector<complex<double> > psi3(p,1);
  Norm(psi3,'y');
  cout<<"# s"<<setw(20)<<"p"<<setw(20)<<"Energy"<<endl;
  for(int i = 0; i <= step; i++) {
    s = (double) i / step;
    cout<<setprecision(4)<<s;
    cout<<setprecision(10);
    productevolution(psi3, dt, s);
    //Magnetization(psi3);
    successprob(psi3, finalgroundstate);
    Energy(psi3,s);
    norm = Norm(psi3,'n');
    cout<<setw(25)<<setprecision(15)<<norm;
    cout<<endl;
  }
  //std::cerr << "steps here2 " << step << endl;
}
