#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<time.h>
#include<algorithm>
#include<vector>
#include "jacobi.h"

using namespace std;
using namespace arma;
ofstream ofile;

//main function 
int main()
{
	int n, interact;
	//input of the value of n
	cout << "Introduce the value of n (nxn matrix): ";
	cin >> n;
    double tol=1.0E-10, pmin=0, pmax=10, wr=0.0; //initialize tolerance, minimum rho, frequency and max rho
	//input of the case to study and the values of frequency and max rho (if they are needed)
	cout << "Introduce the value of interaction (0 for buckling beam, 1 for harmonic oscillator, 2 for two electrons): ";
	cin >> interact;
	if (interact == 0) {
		pmax = 1.0;
	}
	else if (interact == 1) {
		cout << "Introduce the value of pmax: ";
		cin >> pmax;
	}
	else if (interact == 2) {
		cout << "Introduce the value of frequency: ";
		cin >> wr;
		cout << "Introduce the value of pmax: ";
		cin >> pmax;
	}
	double h = (pmax - pmin) / (double(n)); //computation of the step size
	//n = n - 1;
	clock_t start, end;
	//initialization of matrices and row and column values
	int p,q;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
	initialize(n, h, a, r, v, interact, wr);
	int iterations = 0;
	double max=1.0;
	start = clock();
    //do jacobi algorithm until convergence
	while ( max > tol) {
		offdiag(a, p, q, n, max);
		jacobirot(a, v, p, q, n);
		iterations++;
	}
	end = clock();
	//writing of the time and iterations results
	cout << "Diagonalization took " << iterations << " iterations" << endl;
	cout << scientific << "CPU time (sec) : " << ((double)end - (double)start) / CLOCKS_PER_SEC << endl;
    
    //get analytical and numerical eigenvalue vectors
    vector<double>eigen=get_eigenvals(a,n);
	vector<double>exact(n);
	if (interact == 0) {
		double pi = acos(-1.0);
		for (int i = 0; i < n; i++) {
			exact[i] = 2 / (h*h) - 2 / (h*h)*cos((i+1)*pi / (n + 1));
		}
		//writes the results for the eigenvalues in an output file
		string fileout;
		cout << "Please introduce the name of the output file: ";
		cin >> fileout;
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		ofile << "        Approx:         Exact:      Rel. error: " << endl;
		for (int i = 0; i < n; i++) {
			double RelativeError = fabs((exact[i] - eigen[i]) / eigen[i]);
			ofile << setw(15) << setprecision(8) << eigen[i];
			ofile << setw(15) << setprecision(8) << exact[i];
			ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
		}
		ofile.close();
	}
	else if (interact == 1) {
		exact[0] = 3;
		for (int i = 1; i < n; i++) {
			exact[i] = exact[i - 1] + 4;
		}
		//writes the results for the eigenvalues in an output file
		string fileout;
		cout << "Please introduce the name of the output file: ";
		cin >> fileout;
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		ofile << "        Approx:         Exact:      Rel. error: " << endl;
		for (int i = 0; i < n; i++) {
			double RelativeError = fabs((exact[i] - eigen[i]) / eigen[i]);
			ofile << setw(15) << setprecision(8) << eigen[i];
			ofile << setw(15) << setprecision(8) << exact[i];
			ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
		}
		ofile.close();
	}
	else {
		//writes the results for the eigenvalues in an output file
		string fileout;
		cout << "Please introduce the name of the output file: ";
		cin >> fileout;
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		ofile << "     Eigenvalues:" << endl;
		for (int i = 0; i < n; i++) {
			double RelativeError = fabs((exact[i] - eigen[i]) / eigen[i]);
			ofile << setw(15) << setprecision(8) << eigen[i] << endl;
		}
	}
	ofile.close();
    return 0;
}
