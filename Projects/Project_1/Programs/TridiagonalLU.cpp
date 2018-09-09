#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "time.h"
// use namespace for output and input
using namespace std;
using namespace arma;

// object for  output files
ofstream ofile;
// Functions used
inline double f(double x) {
	return 100.0*exp(-10.0*x);
}
inline double exact(double x) { return 1.0 - (1 - exp(-10))*x - exp(-10 * x); }

// Beginning of the main program
int main(int argc, char *argv[]) {
	clock_t start, finish;  //  declare start and final time
	
	int exponent;
	string filename;
	// We read also the basic name for the output file and the highest power of 10^n we want
	if (argc <= 1) {
		cout << "Bad Usage: " << argv[0] <<
			" read also file name on same line and max power 10^n" << endl;
		exit(1);
	}
	else {
		filename = argv[1]; // first command line argument after name of program
		exponent = atoi(argv[2]);
	}
	// Loop over powers of 10
	for (int i = 1; i <= exponent; i++) {
		start = clock();
		int  n = (int)pow(10.0, i);
		// Declare new file name
		string fileout = filename;
		// Convert the power 10^i to a string
		string argument = to_string(i);
		// Final filename as filename-i-
		fileout.append(argument);
		double h = 1.0 / (n);
		double hh = h * h;
		n = n - 1; // This way only points between endpoints are studied
		// Set up matrix A and vectors b and x
		mat A = zeros<mat>(n, n); vec b(n); vec x(n);
		A(0, 0) = 2.0;  A(0, 1) = -1;  x(0) = h;  b(0) = hh * f(x(0));
		x(n - 1) = x(0) + (n - 1)*h; b(n - 1) = hh * f(x(n - 1));
		
		 //Quick setup of matrix and vectors
		for (int i = 1; i < n-1; i++) {
			A(i, i - 1) = -1.0;
			A(i, i) = 2.0;
			A(i, i + 1) = -1.0;
			x(i) = x(i - 1) + h;
			b(i) = hh * f(x(i));
		}
		A(n - 1, n - 1) = 2.0; A(n - 2, n - 1) = -1.0; A(n - 1, n - 2) = -1.0;
		vec solution = solve(A, b); //Solving of the equation system Ax=b
		finish = clock(); // Finishing time for the algorythm
		double time = ((finish - start) / CLOCKS_PER_SEC); // Conversion of time into second
		cout << "The elapsed time for n = " << n+1 << " was " << time << " seconds." << endl; // Output of time for each exponent
		//Writing of the results in file
		ofile.open(fileout);
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		ofile << "       x:             approx:          exact:       relative error" << endl;
		for (int i = 1; i < n; i++) {
			double xval = x(i);
			double RelativeError = fabs((exact(xval) - solution[i]) / exact(xval));
			ofile << setw(15) << setprecision(8) << xval;
			ofile << setw(15) << setprecision(8) << solution(i);
			ofile << setw(15) << setprecision(8) << exact(xval);
			ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
		}
		ofile.close();
	}

	return 0;
}
