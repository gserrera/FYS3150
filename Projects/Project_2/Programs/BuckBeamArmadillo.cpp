//  Diagonalizing tridiagonal Toeplitz matrix  with armadillo
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include<time.h>

using namespace  std;
using namespace  arma;

// Begin of main program   

int main()
{
	int       i, n;
	double    RMin, RMax, Step, Diag, Nondiag;
	RMin = 0.0; RMax = 1.0;
	cout << "Introduce the value of n (nxn matrix): ";
	cin >> n;
	mat Hamiltonian = zeros<mat>(n, n);
	// Integration step length
	Step = RMax / n;
	Diag = 2.0 / (Step*Step);
	Nondiag = -1.0 / (Step*Step);
	clock_t start, end;
	// Setting up tridiagonal matrix and diagonalization using Armadillo
	Hamiltonian(0, 0) = Diag;
	Hamiltonian(0, 1) = Nondiag;
	for (i = 1; i < n - 1; i++) {
		Hamiltonian(i, i - 1) = Nondiag;
		Hamiltonian(i, i) = Diag;
		Hamiltonian(i, i + 1) = Nondiag;
	}
	Hamiltonian(n - 1, n - 2) = Nondiag;
	Hamiltonian(n - 1, n - 1) = Diag;
	// diagonalize and obtain eigenvalues
	vec Eigval(n);
	start = clock();
	eig_sym(Eigval, Hamiltonian);
	end = clock();
	cout << scientific << "CPU time (sec) : " << ((double)end - (double)start) / CLOCKS_PER_SEC << endl;
	double pi = acos(-1.0);
	cout << "RESULTS:" << endl;
	cout << setiosflags(ios::showpoint | ios::uppercase);
	cout << "Number of Eigenvalues = " << setw(15) << n << endl;
	cout << "Numerical:      Exact:         Error:" << endl;
	for (int i = 0; i < n; i++) {
		double Exact = Diag + 2 * Nondiag*cos((i + 1)*pi / (n + 1));
		cout << setw(15) << setprecision(8) << Eigval[i];
		cout << setw(15) << setprecision(8) << Exact;
		cout << setw(15) << setprecision(8) << fabs(Exact - Eigval[i]) << endl;
	}
	int oioio;
	cout << "Press any key and then enter to exit ";
	cin >> oioio;
	return 0;
}  //  end of main function