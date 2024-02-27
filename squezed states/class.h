#pragma once
#include <iostream>
#include <cmath>
#include <armadillo>
using namespace std;
using namespace arma;

double factorial(int n) {
	double r = 1;
	for (int i = 1; i <= n; i++) {
		r *= i;
	}
	return r;
}

class state {
public:
	cx_mat phi;//matrix of fock states (often phi is used for them)
	cx_vec scs;//squeezed coherent state
	int nphi = 50;//number of phi (cols)
	int nx = 100;//nubmer of x points (rows)
	double xmax = 10;//distance from 0 to furthest point (it will be 2 times xmax from -xmax to xmax)
	double xgap = 2 * xmax / (nx - 1);//distance from point to point
	double pi = acos(-1);
	vec x;

	void normalize(int n) {//integration by simpson rule
		double integral = 0;
		double fx;
		for (int i = 0; i < nx - 1; i++) {
			fx = xgap / 6 * (real(phi(i, n) * conj(phi(i, n))) + real(phi(i + 1, n) * conj(phi(i + 1, n))) + 2 * (real(phi(i, n) * conj(phi(i, n))) + real(phi(i + 1, n) * conj(phi(i + 1, n)))));
			integral += fx;
		}
		phi.col(n) /= sqrt(integral);
	}
	void hermitemaker() {
		//setting up hermite polynominals
		phi = cx_mat(nx, nphi, fill::zeros);
		phi.col(0) += 1;
		for (int j = 0; j < nx; j++) {
			phi(j, 1) = 2 * x[j] * phi(j, 0);
		}
		for (int i = 2; i < nphi; i++) {
			for (int j = 0; j < nx; j++) {
				phi(j, i) = 2 * x[j] * phi(j, i - 1) - double(2 * (i - 1)) * phi(j, i - 2);//armadillo has problems with int
			}
		}
	}

	void squeeze(cx_double a, double b, double t) {//really squeezing (b) and shifting (a) depending on parameters, t is time for time evolution
		cx_vec outl,outm;
		scs=cx_vec(nx,fill::zeros);

		for (int k = 0; 2 * k < nphi;k++) {

			outm = cx_vec(nx, fill::zeros);
			for (int m = 0; m <= 2*k;m++) {

				outl = cx_vec(nx, fill::zeros);
				for (int l = 0; 2 * k - m + l < nphi;l++) {
					outl += pow(a,l)/factorial(l)*sqrt(factorial(2*k-m+l))/factorial(2*k-m)*exp(-1i*t*(2*k-m+l+0.5))*phi.col(2*k-m+l);
				}
				outl *= pow(conj(a), m) / factorial(m)/sqrt(factorial(2 * k - m));
				outm += outl;
			}
			outm *= pow(-tanh(b), k) * factorial(2 * k) / (pow(2, k) * factorial(k));
			scs += outm;
		}
		scs *= sqrt(1 / cosh(b)) * exp(-a * conj(a) / double(2));
	}

	state() {
		x = vec(nx);
		for (int i = 0; i < nx; i++)
			x[i] = -xmax + i * xgap;

		hermitemaker();
		for (int i = 0; i < nphi; i++) {
			for (int j = 0; j < nx; j++)
				phi(j, i) *= 1 / sqrt(pow(2, i) * factorial(i)) * pow((1 / pi), (1 / 4)) * exp(-x[j] * x[j] / 2);//1 / sqrt(pow(2, i) * factorial(i))(1 / pi) ^ (1 / 4)
		}
		for (int i = 0; i < nphi; i++)
			normalize(i);
	}

};