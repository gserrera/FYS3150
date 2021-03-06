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

	cout << "Full solar system (Velocity Verlet)" << endl;
	Dimension = 2;

	IntegrationPoints = 30000;
	FinalTime = 300.;

	double TimeStep = FinalTime / ((double)IntegrationPoints); //Calculates the time step
	// initial position x = 1AU, y = z = 0, vx = 2pi, vy=0, vz=0

	//DEFINE ALL PLANETS AND SOLVER SYSTEM

	planet planet1(0.000003, 0.935, 0.354, -2.28, 5.87); // Earth: (relative mass,x,y,vx,vy)
	planet planet2(1., -0.000147, 0.00725, -0.00223, 0.0011);           // Sun
	planet planet3(0.00000033, 1.37, 0.15, 0.77, 5.51); //Mars
	planet planet4(0.00095, -2.66, -4.66, 2.36, -1.23); //Jupiter
	planet planet5(0.000275, 1.55, -9.93, 1.898, 0.31); //Saturn
	planet planet6(0.000000165, -0.15, -0.43, 7.62, -2.87); //Mercury
	planet planet7(0.00000245, -0.7, -0.17, -1.66, 7.17); //Venus
	planet planet8(0.000044, 17.2, 9.99, -0.73, 1.17); //Uranus
	planet planet9(0.0000515, 28.92, -7.72, 0.29, 1.11); //Neptune
	solver solarsystem_vv;

	//Add planets to solver

	solarsystem_vv.add(planet2);
	solarsystem_vv.add(planet6);
	solarsystem_vv.add(planet7);
	solarsystem_vv.add(planet1);
	solarsystem_vv.add(planet3);
	solarsystem_vv.add(planet4);
	solarsystem_vv.add(planet5);
	solarsystem_vv.add(planet8);
	solarsystem_vv.add(planet9);


    //Solving of the ODE's and simulation of the whole system
	solarsystem_vv.VelocityVerlet(Dimension, IntegrationPoints, FinalTime, 9);
	
	
	return 0;
}