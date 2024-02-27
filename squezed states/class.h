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

class squeeze {
public:
	cx_mat phi;
	int nphi = 10;//number of phi (cols)
	int nx = 100;//nubmer of x points (rows)
	double xmax = 10;//distance from 0 to furthest point (it will be 2 times xmax from -xmax to xmax)
	double xgap = 2 * xmax / (nx - 1);//distance from point to point
	double pi = acos(-1);
	vec x;

	squeeze() {
		x = vec(nx);
		for (int i = 0; i < nx; i++)
			x[i] = -xmax + i * xgap;

		//setting up hermite polynominals
		phi = cx_mat(nx, nphi, fill::zeros);
		phi.col(0) += 1;
		for (int j = 0; j < nx; j++) {
			phi(j,1) = 2 * x[j] * phi(j, 0);
		}
		for (int i = 2; i < nphi; i++) {
			for (int j = 0; j < nx; j++){
				phi(j, i) = 2 * x[j] * phi(j, i-1)- double(2 * (i - 1)) * phi(j, i-2);//armadillo has problems with int
			}
		}
		// time for wavefunctions
		for (int i = 0; i < nphi; i++) {
			for (int j = 0; j < nx; j++)
				phi(j,i)*= 1 / sqrt(pow(2, i) * factorial(i))*pow((1 / pi), (1 / 4))*exp(-x[j]*x[j]/2);//1 / sqrt(pow(2, i) * factorial(i))(1 / pi) ^ (1 / 4)
		}
}

};