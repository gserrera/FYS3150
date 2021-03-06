// This code uses standard c++ allocation of arrays

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include "planet.h"
#include "solver.h"


using namespace std;



int main()
{
	int IntegrationPoints;  // Integration points
	double FinalTime;       // End time of calculation (in years)
	int Dimension;           // No. of spatial dimensions

	cout << "Earth-Sun-Jupiter system (Velocity Verlet)" << endl;
	Dimension = 2;
	int n;
	cout << "Introduce factor for Jupiter's mass: ";
	cin >> n;
	IntegrationPoints = 10000;
	FinalTime = 100.;

	double TimeStep = FinalTime / ((double)IntegrationPoints); //Calculates the time step

	//DEFINE ALL PLANETS AND SOLVER SYSTEM

	planet planet1(0.000003, 0.935, 0.347, -2.28, 5.84); // Earth: (relative mass,x,y,vx,vy)
	planet planet2(1., 0.0, 0.0, 0.0, 0.0);           // Sun
	planet planet3(0.00095*n, -2.66, -4.66, 2.36, -1.23); // Jupiter

	
	solver earthsun_vv;

	//Add planets to solver

	earthsun_vv.add(planet2);
	earthsun_vv.add(planet1);
	earthsun_vv.add(planet3);


	
    //Solving of the ODE's and simulation of the whole system
	earthsun_vv.VelocityVerlet(Dimension, IntegrationPoints, FinalTime, 3);
	
	return 0;
}