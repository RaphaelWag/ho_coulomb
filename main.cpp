//Author: Raphael Wagner 20.09.2018

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono> //used to measure runtime
#include "Jacobi_method.h"

void get_lowest_eigenvalue(double **A, double &lowest_eigenvalue, int N);
void fill_matrices(double **A, double **R, double potential(double,double), int N, double h, double omega_2);
double potential(double x, double omega_2);

using namespace std;
using namespace std::chrono;

int main()
{
    int N = 150;
    int max_values = 3;
    double rho_max_array[] = {42.5,45.0,47.5};
    double lowest_eigenvalue = 0;
    int omegavalues = 5;
    double omega = 0;
    double omega_array[]= {0.01,0.25,0.5,1.0,5.0};
    double eigenvalues[3];
    double mean;
    double std;


    for (int ii = 0; ii < omegavalues; ++ii){

        omega = omega_array[ii];
        double omega_2 = omega*omega;
        cout << "w_r:" << omega << endl;

            for (int n = 0; n < max_values; ++n) {

                double rho_max = (rho_max_array[n])*pow(0.8,ii);
                if(ii==omegavalues-1){
                    rho_max = rho_max*0.8;
                }

                cout << "rho_max=" << rho_max << ":" << endl;

                //set calculation parameter

                double h = rho_max / double(N+1); //stepsize

                int rotations = 0;
                double tolerance = pow(10, -16);

                //construct arrays

                auto **R = new double *[N]; //identity matrix to store Eigenvectors
                for (int i = 0; i < N; ++i) {
                    R[i] = new double[N];
                }

                auto **A = new double *[N];
                for (int i = 0; i < N; ++i) {
                    A[i] = new double[N];
                }

                //set array values

                fill_matrices(A,R,potential,N,h,omega_2);

                //Perform Jacobi Rotation

                //start time
                auto start = high_resolution_clock::now();

                Jacobi_EV(tolerance,rotations,A,R,N);

                auto stop = high_resolution_clock::now();
                auto duration = duration_cast<nanoseconds>(stop - start);
                auto time = duration.count();

                get_lowest_eigenvalue(A, lowest_eigenvalue, N);

                //print results

                cout << "lowest eigenvalue:" << lowest_eigenvalue << endl;
                cout << "runtime / seconds:" << time * pow(10, -9) << endl;
                cout << endl;

                eigenvalues[n] = lowest_eigenvalue;

                //destruct arrays

                for (int i = 0; i < N; ++i) {
                    delete[] R[i];

                }
                delete[] R;

                for (int i = 0; i < N; ++i) {
                    delete[] A[i];

                }
                delete[] A;
            }
            mean = eigenvalues[0]+eigenvalues[1]+eigenvalues[2];
            mean = mean/3.0;
            std = sqrt(pow(mean-eigenvalues[0],2)+pow(mean-eigenvalues[1],2)+pow(mean-eigenvalues[2],2))/sqrt(3);
            cout << "statistical lowest eigenvalue:" << mean << "+-" << std << endl;
            cout << endl;

    }
    return 0;
}

void get_lowest_eigenvalue(double **A, double &lowest_eigenvalue, int N)
{
    lowest_eigenvalue = A[0][0];

    for (int i = 0; i < N; ++i)
    {
        if (lowest_eigenvalue > A[i][i])
            lowest_eigenvalue = A[i][i];
    }
}

void fill_matrices(double **A, double **R, double potential(double,double), int N, double h, double omega_2)
{
    double hh = h*h;

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j)
            {
                R[i][j] = 1;
                A[i][j] = 2.0/hh + potential((i+1)*h, omega_2);
            } else
            {
                R[i][j] = 0;
                A[i][j] = 0;
            }
            if ((i == (j + 1)) or (i == (j - 1)))
            {
                A[i][j] = -1.0/hh;
            }
        }
    }
}

double potential(double x, double omega_2)
{
    //harmonic oscillator and coulomb interaction
    return omega_2 * x*x + 1.0/x;
}
