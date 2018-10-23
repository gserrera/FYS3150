#include "catch.hpp"
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

#define PI 3.14159265
TEST_CASE("Testing conservation of energy and angular momentum") {
	int IntegrationPoints;  // Integration points
	double FinalTime;       // End time of calculation (in years)
	int Dimension;           // No. of spatial dimensions
	Dimension = 2;
	FinalTime = 100.;
	for (int n = 0; n <= 6; n++) { //Changes the number of steps and therefore the step size
		IntegrationPoints = pow(10,n);
		double TimeStep = FinalTime / ((double)IntegrationPoints); //Calculates the time step

		//DEFINE ALL PLANETS AND SOLVER SYSTEM

		planet planet1(0.000003, 1.0, 0.0, 0.0, 6.3); // Earth: (relative mass,x,y,vx,vy)
		planet planet2(1., 0.0, 0.0, 0.0, 0.0);           // Sun


		solver earthsun_vv;

		//Add planets to solver

		earthsun_vv.add(planet2);
		earthsun_vv.add(planet1);

		double K1 = planet1.KineticEnergy(); //Compute Energies and momentum at the start
		double V1 = planet1.PotentialEnergy(planet2, 4 * PI*PI);
		double L1 = sqrt(planet1.velocity[0] * planet1.velocity[0] + planet1.velocity[1] * planet1.velocity[1])*sqrt(planet1.position[0] * planet1.position[0] + planet1.position[1] * planet1.position[1]);
		earthsun_vv.VelocityVerlet(Dimension, IntegrationPoints, FinalTime, 2);

		double K2 = planet1.KineticEnergy(); //Compute energies and momentum at the end
		double V2 = planet1.PotentialEnergy(planet2, 4 * PI*PI);
		double L2 = sqrt(planet1.velocity[0] * planet1.velocity[0] + planet1.velocity[1] * planet1.velocity[1])*sqrt(planet1.position[0] * planet1.position[0] + planet1.position[1] * planet1.position[1]);
		
		REQUIRE(K2 == Approx(K1).epsilon(0.00001));
		REQUIRE(V2 == Approx(V1).epsilon(0.00001));
		REQUIRE(L2 == Approx(L1).epsilon(0.00001));
	}
}