#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream ofile;

int main()
{
	string filename;
	cout << "Introduce the name of the output file: ";
	cin >> filename;
	int n = 99;
	int tsteps = 40000;
	double alpha = (n+1)*(n+1) /((double)tsteps*tsteps);
	//Creates matrices for allocating the results (in 2D) in a given moment and the next one
	mat u(n + 1, n + 1); mat unew(n + 1, n + 1);
	for (int i = 0; i <= n; i++) { //Boundary conditions
		u(i, 0) = unew(i, 0) = u(i, n) = unew(i, n) = u(n, i) = unew(n, i) = 0.0;
		u(0, i) = unew(0, i) = 1.0;
	}
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			//  initial condition
			u(i,j) = 0.0;
			//  intitialise the new matrix 
			unew(i,j) = 0.0;
		}
	}
	string fileout = filename;
	ofile.open(fileout);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "       time:           x:          y:         u:" << endl;
	cout << "The value of alpha is " << alpha;
	// Time integration
	for (int t = 1; t <= tsteps; t++) {
		for (int i = 2; i < n; i++) {
			for (int j = 2; j < n; j++) {
				// Discretized diff equation and prints results
				unew(i, j) = (1 - 4 * alpha) * u(i, j) + alpha * (u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1));

				ofile << setw(15) << setprecision(8) << t;
				ofile << setw(15) << setprecision(8) << i;
				ofile << setw(15) << setprecision(8) << j;
				ofile << setw(15) << setprecision(8) << unew(i,j) << endl;
			}
		}
		for (int i2 = 1; i2 < n; i2++) {
			for (int j2 = 1; j2 < n; j2++) {
				u(i2, j2) = unew(i2, j2); //Replaces old matrix with new for the next time step
			}
		}
		
	}
	ofile.close();  // close output file
	return 0;
}
