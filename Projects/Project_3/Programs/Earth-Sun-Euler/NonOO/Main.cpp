#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace  std;

ofstream ofile;

void output(double, double, double, double, double);

int main() //Main function for the Sun-Earth system using a non-object oriented Euler method
{
	cout << "Sun-Earth system (Euler method, not Object oriented)" << endl;
	//Declaration of variables
	string fileout;
	cout << "Please introduce the name of the output file: ";
	cin >> fileout;

	ofile.open(fileout);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "        Time:        x:           y:             vx:              vy: " << endl;

	int NumberofSteps = 1000;
	double FinalTime = 1.0;
	double Step = FinalTime / ((double)NumberofSteps);
	double time = 0.0;
	//Initial values  x = 1.0 AU and vy = 2*pi
	double pi = acos(-1.0);
	double FourPi2 = 4 * pi*pi;
	double x = 1.0; double y = 0.0; double vx = 0.0; double vy = 2.0*pi;
	double r = sqrt(x*x + y*y);

	//Solving of the ODE's using Euler method
	while (time <= FinalTime) {
		x += Step * vx;
		y += Step * vy;
		vx -= Step * FourPi2*x / (r*r*r);
		vy -= Step * FourPi2*y / (r*r*r);
		r = sqrt(x*x + y * y);
		time += Step;
		output(time, x, y, vx, vy);   // Write the results to file 
	}
	ofile.close(); 
	return 0;
}   //End of main function 

//Function to write out the results to the output file
void output(double time, double x, double y, double vx, double vy)
{
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setw(15) << setprecision(8) << time;
	ofile << setw(15) << setprecision(8) << x;
	ofile << setw(15) << setprecision(8) << y;
	ofile << setw(15) << setprecision(8) << vx;
	ofile << setw(15) << setprecision(8) << vy << endl;
}  // end of function output
