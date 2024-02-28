#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <armadillo>//great vectors and complex numbers compared to <vector> and <complex>
#include "gnuplot.h"

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
	int nphi = 80;//number of phi (cols)
	int nx = 200;//nubmer of x points (rows)
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
	void plotstate(cx_double a, double b) {//not sure how to do it without using a file so it's a dumb way
		ofstream file("data-a_" + to_string(float(real(a))) + " i" + to_string(float(imag(a))) + "b_" + to_string(float(b)) + ".txt");
		GnuplotPipe gp;
		string s;
		mat output(nx,3);
		double tmax = 4 * pi;
		int N=100;
		for (int i = 0; i < N; i++) {
			cout << "step " << i + 1 << " out of " << N << endl;
			squeeze(a, b, i*tmax/(N-1));
			output.col(0) = x;
			output.col(1) = real(scs);
			output.col(2) = imag(scs);
			file << output << endl<<endl;
		}
		gp.sendLine("set xrange [-10:10];");
		gp.sendLine("set yrange [-1.5:1.5];");
		gp.sendLine("set terminal gif animate delay 8");
		gp.sendLine("set output \"wavefunction-a_" + to_string(float(real(a))) + " i" + to_string(float(imag(a))) + "b_" + to_string(float(b)) + ".gif\"");
		gp.sendLine("do for [i = 1:" + to_string(N) + "] {");
		gp.sendLine("plot \"data.txt\" using 1:2 index i-1 title 'real' with lines, \"data.txt\" using 1:3 index i-1  title 'imag' with lines");
		gp.sendLine("}");
		file.close();
	}
	void plotstatems(cx_double a, double b) {//not sure how to do it without using a file so it's a dumb way
		ofstream file("data-a_" + to_string(float(real(a))) + " i" + to_string(float(imag(a))) + "b_" + to_string(float(b)) + ".txt");
		GnuplotPipe gp;
		string s;
		mat output(nx, 3);
		double tmax = 4 * pi;
		int N = 100;
		for (int i = 0; i < N; i++) {
			cout << "step " << i + 1 << " out of " << N << endl;
			squeeze(a, b, i * tmax / (N - 1));
			output.col(0) = x;
			output.col(1) = real(conj(scs)%scs);
			file << output << endl << endl;
		}
		gp.sendLine("set xrange [-10:10];");
		gp.sendLine("set yrange [-1.5:1.5];");
		gp.sendLine("set terminal gif animate delay 8");
		gp.sendLine("set output \"modulosquared-a_" + to_string(float(real(a))) + " i" + to_string(float(imag(a))) + "b_" + to_string(float(b)) + ".gif\"");
		gp.sendLine("do for [i = 1:" + to_string(N) + "] {");
		gp.sendLine("plot \"data.txt\" using 1:2 index i-1 title 'mod^2' with lines");
		gp.sendLine("}");
		file.close();
	}

	void squeeze(cx_double a, double b, double t) {//really squeezing (b) and shifting (a) depending on parameters, t is time for time evolution
		cx_vec outl,outm;//temporary
		scs=cx_vec(nx,fill::zeros);

		for (int k = 0; 2 * k < nphi;k++) {

			outm = cx_vec(nx, fill::zeros);
			for (int m = 0; m <= 2*k;m++) {

				outl = cx_vec(nx, fill::zeros);
				for (int l = 0; 2 * k - m + l < nphi;l++) {
					outl += pow(a,l)/factorial(l)*sqrt(factorial(2*k-m+l))*exp(-1i*t*(2*k-m+l+0.5))*phi.col(2*k-m+l);
				}
				outl *= pow(-conj(a), m) / factorial(m)/factorial(2 * k - m);
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