#ifndef PLANET_H
#define PLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;


class planet
{
public:
	// Properties
	double mass;
	double position[2];
	double velocity[2];
	double potential;
	double kinetic;

	// Initializers
	planet();
	planet(double M, double x, double y, double vx, double vy);

	// Functions
	double distance(planet otherPlanet);
	double GravitationalForce(planet otherPlanet, double Gconst);
	double KineticEnergy();
	double PotentialEnergy(planet &otherPlanet, double Gconst);

};

#endif // PLANET_H#pragma once
