#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream ofile;

int main() //Main program. Prints result matrix for two different times.
{
	string filename; string filename2;
	cout << "Introduce the name of the output file for a small time result: ";
	cin >> filename;
	cout << "Introduce the name of the output file for a large time result: ";
	cin >> filename2;
	int n = 99;
	int tsteps = 40000;
	int tsmall = tsteps / 1000;
	int tlarge = tsteps * 3 / 4;
	double alpha = (n+1)*(n+1)/ ((double)tsteps);

	double step = 1 / n;
	mat u(n + 1, n + 1); mat unew(n + 1, n + 1);
	for (int i = 0; i <= n; i++) { // boundary conditions
		u(i, 0) = unew(i, 0) = u(i, n) = unew(i, n) = u(n, i) = unew(n, i) = 0.0;
		u(0, i) = unew(0, i) = 1.0;
	}
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			//  initial condition
			u(i, j) = 0.0;
			//  intitialise the new matrix
			unew(i, j) = 0.0;
		}
	}
	string fileout = filename;
	string fileout2 = filename2;

	cout << "The value of alpha is " << alpha;
	// Time integration
	for (int t = 1; t <= tsteps; t++) {
		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				// Discretized diff eq
				unew(i, j) = (1 - 4 * alpha) * u(i, j) + alpha * (u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1));

			}
		}
		for (int i2 = 1; i2 < n; i2++) { //Copy new matrix into old
			for (int j2 = 1; j2 < n; j2++) {
				u(i2, j2) = unew(i2, j2);
			}
		}
		if (t == tsmall) { //Prints matrix for two different time steps.
			ofile.open(fileout);
			ofile << setiosflags(ios::showpoint | ios::uppercase);
			ofile << " Small time matrix" << endl;
			for (int i3 = 0; i3 <= n; i3++) {
				for (int j3 = 0; j3 <= n; j3++) {
					ofile << setw(15) << setprecision(8) << unew(i3, j3);
				}
				ofile << endl;
			}
			ofile.close();  // close output file
		} else if (t == tlarge) {
			ofile.open(fileout2);
			ofile << setiosflags(ios::showpoint | ios::uppercase);
			ofile << " Large time matrix" << endl;
			for (int i3 = 0; i3 <= n; i3++) {
				for (int j3 = 0; j3 <= n; j3++) {
					ofile << setw(15) << setprecision(8) << unew(i3, j3);
				}
				ofile << endl;
			}
			ofile.close();  // close output file
		}
	}

}
