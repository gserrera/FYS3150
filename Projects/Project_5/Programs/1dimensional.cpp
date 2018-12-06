#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include<algorithm>
#include<vector>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream ofile;

void forward_euler(int, int, double, string);
void backward_euler(int, int, double, string);
void crank_nicolson(int, int, double, string);
void tridiag(double, double, vec&, vec&, int);

//Beginning of main function
int main()
{
	string filename;
	int n; int tsteps;
	cout << "Introduce the name of the output file: ";
	cin >> filename;
	cout << "Introduce the number of x steps: ";
	cin >> n;
	cout << "Introduce the number of time steps: ";
	cin >> tsteps;
	int method;
	cout << "Introduce the solving method (0 for explicit, 1 for implicit, 2 for Crank-Nicolson): ";
	cin >> method;

	if (method == 0) {
		n = n - 1;
		double alpha = (n + 1)*(n + 1) / (double)(tsteps);
		cout << "The value of alpha is " << alpha << endl;
		forward_euler(n, tsteps, alpha, filename);
		cout << "Done";
	}
	else if (method == 1) {
		n = n - 3;
		double alpha = (n + 1)*(n + 1) / (double)(tsteps);
		cout << "The value of alpha is " << alpha << endl;
		backward_euler(n, tsteps, alpha, filename);
		cout << "Done";
	}
	else if (method == 2) {
		n = n - 3;
		double alpha = (n + 1)*(n + 1) / (double)(tsteps);
		cout << "The value of alpha is " << alpha << endl;
		crank_nicolson(n, tsteps, alpha, filename);
		cout << "Done";
	}
	
	return 0;
} 
void forward_euler(int n, int tsteps, double alpha, string filename){

	vec u(n + 1); vec unew(n + 1);
	u(0) = unew(0) = 0.0;
	u(n) = unew(n) = 1.0;
	for (int i = 1; i < n; i++) {
		//  initial condition
		u(i) = 0.0;
		//  intitialise the new vector 
		unew(i) = 0.0;
	}
	string fileout = filename;
	ofile.open(fileout);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "       time:           x:          u:" << endl;
	// Time integration
	for (int t = 1; t <= tsteps; t++) {
		for (int i = 1; i < n; i++) {
			// Discretized diff eq
			unew(i) = alpha * u(i - 1) + (1 - 2 * alpha) * u(i) + alpha * u(i + 1);
		} for (int i = 0; i <= n; i++) { //Print results
			ofile << setw(15) << setprecision(8) << t;
			ofile << setw(15) << setprecision(8) << i;
			ofile << setw(15) << setprecision(8) << unew(i) << endl;
		}
		for (int j = 1; j < n; j++) u(j) = unew(j); //Copy next step into previous step
	
	}
	ofile.close();  // close output file

}

void backward_euler(int n, int tsteps, double alpha, string filename)
{
	double a, b;
	vec y(n+1); // Right side of matrix equation Au=y, the solution at a previous step
	vec sol(n + 3); //Final solution
	vec u(n + 1); // Left side of matrix equation Au=y
	// Initial conditions
	for (int i = 1; i < n; i++) {
		y(i) = 0.0;
	}
	// Boundary conditions
	y(n) = 1.0;
	y(0) = 0.0;
	// Matrix A, only constants
	a = - alpha;
	b = 1.0 + 2.0 * alpha;
	string fileout = filename;
	ofile.open(fileout);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "       time:           x:          u:" << endl;
	// Time iteration
	for (int t = 1; t <= tsteps; t++) {
		//  here we solve the tridiagonal linear set of equations, 
		tridiag(a, b, y, u, n);
		//Creates the solution vector with the boundary conditions
		sol(0) = 0.0;
		sol(n + 2) = 1.0;
		for (int i = 1; i < n + 2; i++) { //Constructs the final solution and prints it
			sol(i) = u(i - 1);
		}
		// replace previous time solution with new
		for (int i = 0; i <= n; i++) {
			y(i) = u(i);
		}
		y(n) = u(n) - a; //Last element is hardcoded from the boundary condition.
		//  print results
		for (int i = 0; i <= n+2; i++) {
			ofile << setw(15) << setprecision(8) << t;
			ofile << setw(15) << setprecision(8) << i;
			ofile << setw(15) << setprecision(8) << sol(i) << endl;
		}
	}   // end time iteration
	ofile.close();
} 
void tridiag(double a, double b, vec& y, vec& u, int n){

	mat A = zeros<mat>(n + 1, n + 1); //Constructs the A matrix
	A(0, 0) = b;  A(0, 1) = a;
	for (int i = 1; i < n; i++) {
		A(i, i - 1) = a;
		A(i, i) = b;
		A(i, i + 1) = a;
	}
	A(n, n) = b; A(n - 1, n) = a; A(n, n - 1) = a;
	u = solve(A, y); //Solves the set of equations using LU decomposition
	
}
void crank_nicolson(int n, int tsteps, double alpha, string filename)
{
	double a, b;
	vec r(n + 1); // Right side of matrix equation Au=r, the solution at a previous step
	vec u(n + 1); // Left side of matrix equation Au=r
	vec sol(n + 3); //Final solution
	// Initial conditions
	for (int i = 1; i < n; i++) {
		sol(i) = 0.0;
	}
	// Boundary conditions
	sol(n+2) = 1.0;
	sol(0) = 0.0;
	// setting up the matrix 
	a = -alpha;
	b = 2 + 2 * alpha;
	string fileout = filename;
	ofile.open(fileout);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "       time:           x:          u:" << endl;
	for (int i = 0; i <= n; i++) {
		r(i) = alpha * sol(i) + (2 - 2 * alpha)*sol(i+1) + alpha * sol(i + 2);  //Calculate r following the explicit scheme
	}
	// Time iteration
	for (int t = 1; t <= tsteps; t++) {
		tridiag(a, b, r, u, n);
		sol(0) = 0.0;
		sol(n+2) = 1.0;
		for (int i = 1; i < n + 2; i++) { //Constructs the final solution
			sol(i) = u(i - 1);
		}
		for (int i = 0; i <= n+2; i++) { //Printing
			ofile << setw(15) << setprecision(8) << t;
			ofile << setw(15) << setprecision(8) << i;
			ofile << setw(15) << setprecision(8) << sol(i) << endl;
		}
		for (int i = 0; i <=n; i++) {
			r(i) = alpha * sol(i) + (2 - 2 * alpha)*sol(i + 1) + alpha * sol(i + 2); //Calculate r following the explicit scheme
		}
		r(n) += -a; //Hardcode last element with boundary condition
	}
	ofile.close();
}
