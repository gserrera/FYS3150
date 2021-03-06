#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace  std;

ofstream ofile;

void output(double, double, double, double, double);

int main() //Main function for the Sun-Earth system using a non-object oriented Velocity Verlet method
{
	cout << "Sun-Earth system (Velocity Verlet, not Object oriented)" << endl;
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
	double r = sqrt(x*x + y * y);

	//Solving of the ODE's using Velocity Verlet method
	while (time <= FinalTime) {
		double ax1 = -4 * pi*pi*x / (r*r*r);
		double ay1 = -4 * pi*pi*y / (r*r*r);
		x += Step * vx + 0.5*Step*Step*ax1;
		y += Step * vy + 0.5*Step*Step*ay1;
		double ax2 = -4 * pi*pi*x / (r*r*r);
		double ay2 = -4 * pi*pi*y / (r*r*r);
		vx += 0.5*Step*ax1 + 0.5*Step*ax2;
		vy += 0.5*Step*ay1 + 0.5*Step*ay2;
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